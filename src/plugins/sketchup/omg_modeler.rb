##
# OMG, Optimization-based Model Generator
# Copyright (c) 2016-2017 The University of Hong Kong
# Author: Fan Xue <xuef@hku.hk; fanxue@outlook.com>
#
# This file is part of omg.
#
# OMG is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OMG is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with OMG.  If not, see <http://www.gnu.org/licenses/>.
##

begin
  $LOAD_PATH << File.dirname(__FILE__) + "/omg_modeler"
  $LOAD_PATH.uniq!
  Sketchup.send_action('showRubyPanel:')
end

require 'sketchup.rb'
require 'omg'

SKETCHUP_CONSOLE.show

$omg_ground = nil
$omg_temp_folder = "Z:/"

class Omg_model
  def initialize(solvername = nil)
    @solver = Omg::Omg.new()
    @iterations_incremental = 10
    @web_dialog = nil
    @best_fitness = 1000000.0
    @iteration_overall = 0;
    @best_so_far = 1.0E8;
    @best_incremental = 1.0;
    @progress_js_file = $omg_temp_folder + "progress.js"
    @log_file = $omg_temp_folder + "log.txt"
    @temp_image_file = $omg_temp_folder + "1.bmp"
    @best_image_file = $omg_temp_folder + "best.bmp"
    @VERBOSE = false
    @history_cache = 10
    @history = Array.new(@history_cache, 0.0)
    @historyB = Array.new(@history_cache, 0.0)
    @historyI = Array.new(@history_cache, 0.0)
    @message = ""
	@refine_mode = false
	
    open(@log_file, 'w') do |f|
	  f << "\n"
	end
    
    @time_overall = 0.0
    @time_client = 0.0
    @time_evaluation = 0.0
    @time_algorithm = 0.0
    reset_guide_planes()
    # no edge display
    model = Sketchup.active_model
    model.rendering_options["DrawSilhouettes"] = false
    model.rendering_options["ExtendLines"] = false
    model.rendering_options["EdgeDisplayMode"] = 0
    model.active_view.refresh
    UI.refresh_inspectors()
  end
  def set_current_components(comp)
    @current_component = comp
  end
  def get_solver()
    return @solver
  end
  def get_temp_file()
    return @temp_image_file
  end
  def get_best_file()
    return @best_image_file
  end
  def get_verbose()
    return @VERBOSE
  end
  def set_iterations_incremental(iter)
    @iterations_incremental = iter
  end
  def set_ref_measurement(filename)
    # default: SSIM
	# 1: MSSSIM
	# 2: PSNR
    # 4: MSE
    @solver.load_ref_measurement(filename, -9999)
  end
  def set_web_dialog(web_dialog)
    @web_dialog = web_dialog
  end
  def set_geolocation(lat, lng, country = nil, address = nil)
    shadowinfo = Sketchup.active_model.shadow_info
    shadowinfo["Latitude"] = lat
    shadowinfo["Longitude"] = lng
    shadowinfo["Country"] = country || "HK"
    shadowinfo["Location"] = address || "HK"
  end
  def set_sunrise_time(utc_year, utc_month, utc_day, utc_hour)
    shadowinfo = Sketchup.active_model.shadow_info
    firstDay = Date.new(2016,4,5)
    timeZoneOffsetInSecond = firstDay.to_time.utc_offset
    shadowTimeForSunrise = Time.new(utc_year, utc_month, utc_day, utc_hour, 0,0) + timeZoneOffsetInSecond 
    shadowinfo["ShadowTime"] = shadowTimeForSunrise
  end
  def add_guide_planes(array)
    return if array == nil
    array.each do |e| 
      @guide_planes.push(e)
    end
  end
  def remove_guide_plane(idx)
    @guide_planes.delete_at(idx) if idx < @guide_planes.length
  end
  def get_guide_planes()
    return @guide_planes
  end
  def get_iterations_incremental()
    return @iterations_incremental
  end
  def reset_guide_planes()
    @guide_planes = Array.new
  end
  def get_web_dialog()
    return @web_dialog
  end
  def set_iteration_info(value)
    @best_so_far = value if @best_so_far > value
	if @best_incremental > value && !@refine_mode && !@current_component.nil?
      @current_component.record_best_param()
	end
    @best_incremental = value if @best_incremental > value
    @history[@iteration_overall % @history_cache] = value
    @historyI[@iteration_overall % @history_cache] = @best_incremental
    @historyB[@iteration_overall % @history_cache] = @best_so_far
    # log to file
    open(@log_file, 'a') do |f|
      f << @iteration_overall << "," << value << "," << @best_incremental << "," << @best_so_far << "," << @time_overall << "," << @time_client << "," << @time_evaluation << "," << @time_algorithm << "\n";
    end
    @iteration_overall += 1
    if @iteration_overall % @history_cache == 0
      open(@progress_js_file, 'w') do |f|
        f << "var currnet_iteration = " << @iteration_overall << ";\n"
        f << "var omg_message = '" << @message << "';\n"
        f << "var progress = [\n";
        (0..(@history_cache - 1)).each do |idx|
          f << "{'iteration' : " << (@iteration_overall - @history_cache + 1 + idx) << ",\n";
          f << "'best_fitness' : " << 1.0-@historyB[idx] << ",\n";
          f << "'best_incremental' : " << 1.0-@historyI[idx] << ",\n";
          f << "'current_fitness' : " << 1.0-@history[idx] << "},\n";
        end
        f << "];";
      end
    end
  end

  def build_incremental(folder, files, relation, base_obj, manifolds = nil, force_using_guides = nil)
    # equal trials of incremental
    @best_incremental = 1.0
    manifolds ||= 1
    force_using_guides ||= false
	collection = Array.new
    files.each do |file|
      filename = file
      foldername = folder
      if (filename.index('/') != nil)
        foldername = foldername + filename[0, filename.index('/')]
        filename = filename[filename.index('/')+1, filename.length]
      end
      #puts foldername, filename
      comps = Omg_components.new
      vert = Omg_component.new(filename, foldername, base_obj, relation, self, manifolds, force_using_guides)
      comps.add_instance(vert)
	  collection << comps
	  sleep(0.01)
	  comps.hide_all()
    end
	
	bestf = 1E+6
	best_comps = nil
	
    collection.each do |comps|
	  comps.show_all() 
      adapter = Omg_adapter.new(comps, self)
	  comps.set_adapter(adapter)
	  adapter.project_to_file()
      f_before = adapter.evaluate_fitness()
	  #puts "init = #{f_before}"
	  trial_iterations = (@iterations_incremental * manifolds).to_i
      f = adapter.evolve(trial_iterations)
	  adapter.place_component_from_solution(false)
	  if f <= bestf #&& (f <= f_before || manifolds > 1)
	    bestf = f
		best_comps = comps
      end
	  comps.hide_all()
	end
	if best_comps != nil
	  best_comps.show_all()
	  best_comps.get_adapter().place_component_from_best_solution()
      obj_name = "nil"
      begin
        obj_name = base_obj.get_last_component().name()
      rescue
        obj_name = "nil"
      end
      @message = "A component from [#{files.join(',')}] was attached to #{obj_name}."
	end
    return best_comps
  end
  
  def build_refine(parts, enable_two_swap = false)
    parts.each do |obj|
	  obj.reload_num_variables()
      adapter = Omg_adapter.new(obj, self)
	  adapter.project_to_file()
      bestf = adapter.evaluate_fitness()
      f = adapter.evolve(@iterations_incremental)
	  if f <= bestf
	    adapter.place_component_from_solution(true)
      else
	    obj.rollback()
      end
	  adapter.project_to_file()
      nf = adapter.evaluate_fitness()
	  puts "#{bestf} > #{f} = #{nf}"
	end
  end
  def set_refine_mode(ref)
    @refine_mode = ref
  end
  def get_refine_mode()
    return @refine_mode
  end
  
  def add_time(o_time, c_time, e_time, a_time)
    @time_overall += o_time
    @time_client += c_time
    @time_evaluation += e_time
    @time_algorithm += a_time
  end
  
  def redraw_texture()
    model = Sketchup.active_model
    model.rendering_options["Texture"] = false
    model.active_view.refresh
    UI.refresh_inspectors()
    UI.start_timer(1.0, false) {
      Sketchup.active_model.rendering_options["Texture"] = true
      Sketchup.active_model.active_view.refresh
      UI.refresh_inspectors()
    }
    time_rem = @time_overall - @time_client
    c_pct = (@time_client*100/@time_overall).round(2)
    e_pct = (@time_evaluation*100/@time_overall).round(2)
    a_pct = (@time_algorithm*100/@time_overall).round(2)
    rem_pct = (100.0 - c_pct - e_pct - a_pct).round(2)
    puts ">>>> Modeling process finished, fitness=#{@best_so_far.round(4)}.\n\nOverall time: #{@time_overall}s / #{@iteration_overall} iterations\n (SketchUp #{c_pct}%, Evaluation #{e_pct}%, Algorithm: #{a_pct}%, system/disk: #{rem_pct}%)"
  end
end

class Omg_ground
  def initialize(x_min, x_max, y_min, y_max, z)
    @x_min = x_min 
    @x_max = x_max
    @y_min = y_min
    @y_max = y_max
    @z = z
  end
  def bounds(forHorizontal = false)
    boundingbox = Geom::BoundingBox.new
    boundingbox.add([@x_min, @y_min, @z], [@x_max, @y_max, @z])
    return boundingbox
  end
  def get_last_component()
    return self
  end
  def set_freedom(var_name, is_free)
    # skip
  end
  def freedom_value(var_name)
    return 0
  end
  def name()
    return "Ground"
  end
end

class Omg_guide_plane
  def initialize(a, b, c, d)
    # ax +by +cz +d = 0
    @a = a
    @b = b
    @c = c
    @d = d
  end
  
  def vector_to_point(x, y, z)
    # TODO : complete dist
    if @b == 0 && @c == 0
      return (x + @d/@a), 0.0, 0.0
    end
    return 0.0, 0.0, 0.0
  end
end

class Omg_component
  def initialize(filename = nil, folder = nil, base_obj = nil, relation = "ABOVE", omg_model = nil, manifolds = nil, force_using_guides = nil)
    @su_model = Sketchup.active_model
    @instance = nil
    @instances = Array.new
    @base = base_obj.get_last_component()
    @relation = relation
    @model = omg_model
    @manifolds = (manifolds == nil ? 1 : manifolds.to_i)
    @RATIO_TO_GUIDES = 0.5
    @RATIO_TO_GUIDES = 1.000001 if force_using_guides == true
    
    return if filename == nil or folder == nil or @base == nil
      
    # open component
    family_path = Sketchup.find_support_file filename, folder
    family_component = @su_model.definitions.load family_path
    transform0 = (Geom::Transformation.new Geom::Point3d.new 0,0,0)* Geom::Transformation.scaling(1,1,1)
    overall_instance = @su_model.active_entities.add_instance family_component, transform0
    # read the true object
    overall_array = overall_instance.explode
    if overall_array
      overall_array.each do |inst|
        @instance = inst if (inst.get_attribute 'dynamic_attributes', "omg_scale_x_max").to_f > 0
      end
    end
    # if not found
    if @instance == nil
      UI.messagebox "The file '#{filename}' seems not a valid OMG component."
      return
    end
    
    # get scaling range from component
    @instance.move! transform0
	read_semantics()
    
    @instances.push(@instance)
    for i in 2..@manifolds
      new_inst = nil
      overall_instance = @su_model.active_entities.add_instance family_component, transform0
	#puts "#{(overall_instance.get_attribute 'dynamic_attributes', "omg_scale_x_max").to_f}"
      # read the true object
      overall_array = overall_instance.explode
      if overall_array
        overall_array.each do |inst|
	    #puts ">#{(inst.get_attribute 'dynamic_attributes', "omg_scale_x_max").to_f}"
		  # compatible with x-array manifolds
          new_inst = inst if (inst.get_attribute 'dynamic_attributes', "omg_scale_x_max").to_f > 0
        end
      end
      @instances.push(new_inst) if new_inst != nil
    end
    @free_vars = Array.new(7, true)
	init_freedom_check()
  end
  def get_freedom(var_name)
    if var_name == 'scale_x'
	  return @sem_scale_x_max > @sem_scale_x_min && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	elsif var_name == 'scale_y'
	  return @sem_scale_y_max > @sem_scale_y_min && !@sem_lock_aspect_ratio && isVerticalRelation(@relation) && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	elsif var_name == 'scale_z'
	  return @sem_scale_z_max > @sem_scale_z_min && !@sem_lock_aspect_ratio && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	end
    bbox = @base.bounds(!isVerticalRelation(@relation))
	if var_name == 'position_x'
	  return bbox.width > 0 && (@manifolds > 1 || @sem_ignore_guide_planes || @model.get_guide_planes().length < 1) && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	elsif var_name == 'position_y'
	  return (bbox.max.y-bbox.min.y) > 0 && isVerticalRelation(@relation) && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	elsif var_name == 'position_z'
	  return (bbox.max.z-bbox.min.z) > 0 && !isVerticalRelation(@relation) && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	elsif var_name == 'rotation_z'
	  return (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
	end
	# default : scale_x
	return @sem_scale_x_max > @sem_scale_x_min && (@free_vars[freedom_to_id(var_name)] || !@model.get_refine_mode)
  end
  def set_freedom(var_name, is_free)
    @free_vars[freedom_to_id(var_name)] = is_free
  end
  def freedom_value(var_name)
    return (get_freedom(var_name) ? 1 : 0)
  end
  def freedom_to_id(var_name)
    if var_name == 'scale_x'
	  return 0
	elsif var_name == 'scale_y'
	  return 1
	elsif var_name == 'scale_z'
	  return 2
	elsif var_name == 'position_x'
	  return 3
	elsif var_name == 'position_y'
	  return 4
	elsif var_name == 'position_z'
	  return 5
	elsif var_name == 'rotation_z'
	  return 6
	end
	return 0
  end
  def init_freedom_check()
    if @relation == 'ABOVE'
	  set_freedom('position_z', false)
	  @base.set_freedom('scale_z', false)
	  if @sem_equal_width_parent
	    set_freedom('position_y', false)
	    @base.set_freedom('position_y', false)
	  end
	elsif @relation == 'CONTAINS_ON'
	  set_freedom('position_y', false)
	  set_freedom('scale_y', false)
	  @base.set_freedom('position_y', false)
	  @base.set_freedom('scale_y', false)
	end
  end
  def get_manifolds()
    return @manifolds
  end
  def get_model()
    return @model
  end
  def get_instance()
    return @instances[0]
  end
  def get_num_variables()
	
	# scale x, y, z
    @degree_of_freedom = freedom_value('scale_x') + freedom_value('scale_y') + freedom_value('scale_z')
	# position / offset x (span x if manifold), y, z
    @degree_of_freedom += freedom_value('position_x') + freedom_value('position_y') + freedom_value('position_z')
	@degree_of_freedom += freedom_value('rotation_z')
    # manifolds using last freedom
    last_var_using_guides = (@model.get_guide_planes().length > 0 && @manifolds <= 1 && !@sem_ignore_guide_planes ? 1 : 0)
    return @degree_of_freedom + (@manifolds > 1 ? 1 : 0) + last_var_using_guides
  end
  def read_semantics()
    @len_x_org, @len_y_org, @len_z_org = @instance.unscaled_size()
    @sem_scale_x_max = (@instance.get_attribute 'dynamic_attributes', "omg_scale_x_max").to_f
    @sem_scale_x_min = (@instance.get_attribute 'dynamic_attributes', "omg_scale_x_min").to_f
    @sem_scale_y_max = (@instance.get_attribute 'dynamic_attributes', "omg_scale_y_max").to_f
    @sem_scale_y_min = (@instance.get_attribute 'dynamic_attributes', "omg_scale_y_min").to_f
    @sem_scale_z_max = (@instance.get_attribute 'dynamic_attributes', "omg_scale_z_max").to_f
    @sem_scale_z_min = (@instance.get_attribute 'dynamic_attributes', "omg_scale_z_min").to_f
    @sem_rotation_z_max = (@instance.get_attribute 'dynamic_attributes', "omg_rotation_z_max").to_f
    @sem_rotation_z_min = (@instance.get_attribute 'dynamic_attributes', "omg_rotation_z_min").to_f
    @offset_x = (@instance.get_attribute 'dynamic_attributes', "omg_offset_x").to_f
    @offset_y = (@instance.get_attribute 'dynamic_attributes', "omg_offset_y").to_f
    @offset_z = (@instance.get_attribute 'dynamic_attributes', "omg_offset_z").to_f
    @margin_z = (@instance.get_attribute 'dynamic_attributes', "omg_margin_z").to_f
    @sem_ignore_guide_planes = ((@instance.get_attribute 'dynamic_attributes', "omg_ignore_guide_planes") == "true" ? true : false)
    @sem_lock_aspect_ratio = ((@instance.get_attribute 'dynamic_attributes', "omg_lock_aspect_ratio") == "true" ? true : false)
	@sem_equal_width_parent = ((@instance.get_attribute 'dynamic_attributes', "omg_equal_width_parent") == "true" ? true : false)
	
    
    # fit to boundingbox
    bbox = @base.bounds(!isVerticalRelation(@relation))
    if @sem_scale_x_max > @sem_scale_x_min && @len_x_org > 0
      @sem_scale_x_max = [@sem_scale_x_max, bbox.width.to_f / (@len_x_org * @manifolds)].min
      @sem_scale_x_min = [@sem_scale_x_min, bbox.width.to_f / (@len_x_org * @manifolds)].min
    end
    if @sem_scale_y_max > @sem_scale_y_min && @len_y_org > 0 && isVerticalRelation(@relation)
      @sem_scale_y_max = [@sem_scale_y_max, (bbox.max.y-bbox.min.y).to_f / @len_y_org].min
      @sem_scale_y_min = [@sem_scale_y_min, (bbox.max.y-bbox.min.y).to_f / @len_y_org].min
    end
    if @sem_scale_z_max > @sem_scale_z_min && @len_z_org > 0 && !isVerticalRelation(@relation)
      @sem_scale_z_max = [@sem_scale_z_max, (bbox.max.z-bbox.min.z).to_f * (1.0 - @margin_z * 2) / @len_z_org].min
      @sem_scale_z_min = [@sem_scale_z_min, (bbox.max.z-bbox.min.z).to_f * (1.0 - @margin_z * 2) / @len_z_org].min
    end
  end
  def get_guide_planes()
    gplanes = Array.new
    idx = 0
    #if @manifolds <= 1
    #  return nil
    #end
    while idx < @instances.length do
      gplanes.push(Omg_guide_plane.new(1.0, 0.0, 0.0, -@instances.at(idx).bounds.center.x))
      idx +=1
    end
    return gplanes
  end
  def get_force_using_guides()
    return @RATIO_TO_GUIDES >= 1.0
  end
  def record_best_param()
    @gmx, @gmy, @gmz, @gsx, @gsy, @gsz, @grz = @lmx, @lmy, @lmz, @lsx, @lsy, @lsz, @lrz
  end
  def transform(x)
    # move_x, _y, _z; scale_x, _y, _z; rot_z
    # when @manifolds > 1, the first return is padding_x
    mx, my, mz, sx, sy, sz, rz = decode(x)
	@lmx, @lmy, @lmz, @lsx, @lsy, @lsz, @lrz = decode(x)
	if get_model().get_refine_mode()
      mx, my, mz, sx, sy, sz, rz = decode_refine(x)
	end
    if (@manifolds <= 1)
	  # puts "#{@sem_equal_width_parent} #{my} #{sx * @len_x_org} "
	  if rz > 0
        @instance.move! (Geom::Transformation.new Geom::Point3d.new mx, my, mz) * Geom::Transformation.scaling(sx, sy, sz)
		@instance.transform! (Geom::Transformation.rotation(@instance.bounds.center, Z_AXIS, rz.degrees))
	  else
        @instance.move! (Geom::Transformation.new Geom::Point3d.new mx, my, mz) * Geom::Transformation.scaling(sx, sy, sz)
	  end
      if (@sem_equal_width_parent)
        @base.transform_width_with_child(mx, sx * @len_x_org)
      end
    else
	  # manifolds :: non-rotatable
      # base mx last variables
      #UI::messagebox "x#{x[x.length-1]} lenx#{@len_x_org}   padding#{mx}"
      bbox = @base.bounds(!isVerticalRelation(@relation))
      first_mx = x[x.length-1] * (bbox.width.to_f - sx * @len_x_org * @manifolds - (@manifolds - 1) * mx) + bbox.min.x
      #puts "#{first_mx} (#{mx}) #{my} #{mz} #{sx} #{sy} #{sz} / #{x[x.length-1]}"
      delta_x = sx * @len_x_org + mx
      idx = 0
      while idx < @instances.length do
        inst_mx = delta_x * idx + first_mx
        @instances.at(idx).move! (Geom::Transformation.new Geom::Point3d.new inst_mx, my, mz) * Geom::Transformation.scaling(sx, sy, sz)
        idx +=1
      end
    end
  end
  def record_params()
    mybox = @instance.bounds()
    @best_mx = mybox.min().x
    @best_my = mybox.min().y
    @best_mz = mybox.min().z
    @best_sx = (mybox.max().x - mybox.min().x) / @len_x_org
    @best_sy = (mybox.max().y - mybox.min().y) / @len_y_org
    @best_sz = (mybox.max().z - mybox.min().z) / @len_z_org
	@best_sx = 1.0 if @best_sx.nan?
	@best_sy = 1.0 if @best_sy.nan?
	@best_sz = 1.0 if @best_sz.nan?
	@best_transform = @instance.transformation
  end
  def transform_width_with_child(x_child, len_child)
    @instance.move! (Geom::Transformation.new Geom::Point3d.new x_child, @best_my, @best_mz) * Geom::Transformation.scaling(len_child/@len_x_org, @best_sy, @best_sz)
  end
  def name()
    return @instances.name() + (@manifolds > 1 ? " x#{@manifolds}" : "")
  end
  def rollback()
    if (@manifolds <= 1)
	  if @grz > 0
        @instance.move! (Geom::Transformation.new Geom::Point3d.new @gmx, @gmy, @gmz) * Geom::Transformation.scaling(@gsx, @gsy, @gsz)
		@instance.transform! (Geom::Transformation.rotation(@instance.bounds.center, Z_AXIS, @grz.degrees))
	  else
        @instance.move! (Geom::Transformation.new Geom::Point3d.new @gmx, @gmy, @gmz) * Geom::Transformation.scaling(@gsx, @gsy, @gsz)
	  end
      if (@sem_equal_width_parent)
        @base.transform_width_with_child(@gmx, @gsx * @len_x_org)
      end
    else
	#
	end
  end
  def decode(x)
    # reqires @degree_of_freedom float values [0, 1], 
	#puts "#{x[0]} #{x[1]} #{x[2]} >>> "
    for i in 0..(x.length-1)
	  x[i] = 0.0 if x[i] < 0.0
	  x[i] = 1.0 if x[i] > 1.0
    end
    idx = 0
    sx = (get_freedom('scale_x') || (@best_sx.nil?) ? @sem_scale_x_min : @best_sx)
    sy = (get_freedom('scale_y') || (@best_sy.nil?) ? @sem_scale_y_min : @best_sy)
    sz = (get_freedom('scale_z') || (@best_sz.nil?) ? @sem_scale_z_min : @best_sz)
    bbox = @base.bounds(!isVerticalRelation(@relation))
    
    # freedom 1 : scale_x
    if get_freedom('scale_x') && idx < @degree_of_freedom
      # scalable in deginition
      smax = @sem_scale_x_max
      smin = @sem_scale_x_min
      if bbox.width.to_f / (@len_x_org * @manifolds) < smax
        smax = @sem_scale_x_max
      end
      if smin >= smax
        smin = smax
      end if
      sx = x[idx].abs*(smax - smin) + smin
	  #puts "sx = #{sx} "
      idx += 1
    end
    # freedom 2 : scale_y
    if get_freedom('scale_y') && idx < @degree_of_freedom 
      sy = x[idx].abs*(@sem_scale_y_max - @sem_scale_y_min) + @sem_scale_y_min
	  #puts "sy = #{sy} "
      idx += 1
    elsif @sem_lock_aspect_ratio && isVerticalRelation(@relation)
      sy = x[idx - 1].abs*(@sem_scale_y_max - @sem_scale_y_min) + @sem_scale_y_min
    end
    # freedom 3 : scale_z
    if get_freedom('scale_z') && idx < @degree_of_freedom
      sz = x[idx].abs*(@sem_scale_z_max - @sem_scale_z_min) + @sem_scale_z_min
      idx += 1
	  #puts "sz = #{sz} "
    elsif @sem_lock_aspect_ratio
      sz = x[idx - 1].abs*(@sem_scale_z_max - @sem_scale_z_min) + @sem_scale_z_min
    end
    
    # freedom 4 : mx: location of x , or padding on x (manifold)
    mx = (get_freedom('position_x') || (@best_mx.nil?) ? bbox.min.x : @best_mx)
    my = (get_freedom('position_y') || (@best_my.nil?) ? bbox.max.y : @best_my)
    mz = (get_freedom('position_z') || (@best_mz.nil?) ? (isVerticalRelation(@relation) ? bbox.max.z : bbox.min.z) + (bbox.max.z-bbox.min.z) * @margin_z : @best_mz)
    if bbox.width > 0 && idx < @degree_of_freedom && (@manifolds > 1 || @sem_ignore_guide_planes || @model.get_guide_planes().length < 1) && (@free_vars[3] || @model.get_refine_mode())
      if @manifolds <= 1
        # move on x for single
        mx = x[idx] * (bbox.width.to_f - sx * @len_x_org) + bbox.min.x
	    #puts "mx = #{mx} "
      else
        # padding on x for manifolds
        mx = x[idx] * (bbox.width.to_f / (@manifolds - 1.0) - sx * @len_x_org)
      end
      idx += 1
    end
	# freedom 5 : 
    if get_freedom('position_y') && idx < @degree_of_freedom
      my = x[idx] * (bbox.max.y-bbox.min.y).to_f + bbox.min.y
	  #puts "my = #{my} "
      idx += 1
    end
	# freedom 6 : 
    if get_freedom('position_z') && idx < @degree_of_freedom
      mz = [x[idx] * ((bbox.max.z-bbox.min.z) * (1.0 - @margin_z * 2) - sz * @len_z_org), 0].max + bbox.min.z + (bbox.max.z-bbox.min.z) * @margin_z
      idx += 1
    end
	#freedom 7 : rotation z
	rz = 0
    if get_freedom('rotation_z') && idx < @degree_of_freedom
      rz = x[idx] * (@sem_rotation_z_max - @sem_rotation_z_min)
      idx += 1
    end
	
	# freedom 4: case 2
    # lsat variable of snap_to_guide_planes if necessary, with a 50% possibility to resist snapping
    if @model.get_guide_planes.length > 0 && @manifolds <= 1 && !@sem_ignore_guide_planes && x[x.length - 1] < @RATIO_TO_GUIDES
      gplanes = @model.get_guide_planes
      if @manifolds <= 1
        # when not manifolds
        gplane_idx = ([x[x.length - 1].abs / @RATIO_TO_GUIDES, 0.99999999].min * gplanes.length).floor
        # estimated center
        c_x = mx + @offset_x + sx * @len_x_org / 2
        c_y = my + @offset_y + sy * @len_y_org / 2
        c_z = mz + @offset_z + sz * @len_z_org / 2
        # snap to idx
        snap_x, snap_y, snap_z = gplanes.at(gplane_idx).vector_to_point(c_x, c_y, c_z)
        mx -= snap_x
        my -= snap_y
        mz -= snap_z
      else
        # Nothing
      end
    end
    return mx + @offset_x, my + @offset_y, mz + @offset_z, sx, sy, sz, rz
  end
  def decode_refine(x)
    # reqires @degree_of_freedom float values [0, 1], 
    for i in 0..(x.length-1)
	  x[i] = 0.0 if x[i] < 0.0
	  x[i] = 1.0 if x[i] > 1.0
    end
    idx = 0
    sx = @best_sx
    sy = @best_sy
    sz = @best_sz
    bbox = @base.bounds(!isVerticalRelation(@relation))
    
    # freedom 1 : scale_x
    if get_freedom('scale_x') && idx < @degree_of_freedom
      # scalable in deginition
      smax = @best_sx * 1.2
      smin = @best_sx * 0.8
      if bbox.width.to_f / (@len_x_org * @manifolds) < smax
        smax = @sem_scale_x_max
      end
      sx = x[idx].abs*(smax - smin) + smin
      idx += 1
    end
    # freedom 2 : scale_y
    if get_freedom('scale_y') && idx < @degree_of_freedom 
      sy = (x[idx].abs * 0.4 + 0.8) * @best_sy
	  #puts "sy = #{sy} "
      idx += 1
    elsif @sem_lock_aspect_ratio && isVerticalRelation(@relation)
      sy = (x[idx - 1].abs * 0.4 + 0.8) * @best_sy
    end
    # freedom 3 : scale_z
    if get_freedom('scale_z') && idx < @degree_of_freedom
      sy = (x[idx].abs * 0.4 + 0.8) * @best_sz
      idx += 1
	  #puts "sz = #{sz} "
    elsif @sem_lock_aspect_ratio
      sy = (x[idx - 1].abs * 0.4 + 0.8) * @best_sz
    end
    
    # freedom 4 : mx: location of x , or padding on x (manifold)
    mx = @best_mx - @offset_x
    my = @best_my - @offset_y
    mz = @best_mz - @offset_z
    if bbox.width > 0 && idx < @degree_of_freedom && (@manifolds > 1 || @sem_ignore_guide_planes || @model.get_guide_planes().length < 1) && (@free_vars[3] || @model.get_refine_mode())
      if @manifolds <= 1
        # move on x for single
        mx = (x[idx] * 80 - 40) + @best_mx
	    #puts "mx = #{mx} "
      else
        # padding on x for manifolds
      end
      idx += 1
    end
	# freedom 5 : 
    if get_freedom('position_y') && idx < @degree_of_freedom
      my = (x[idx] * 80 - 40) + @best_my
	  #puts "my = #{my} "
      idx += 1
    end
	# freedom 6 : 
    if get_freedom('position_z') && idx < @degree_of_freedom
      mz = (x[idx] * 80 - 40) + @best_mz
      idx += 1
    end
	#freedom 7 : rotation z
	rz = 0
    if get_freedom('rotation_z') && idx < @degree_of_freedom
      rz = x[idx] * (@sem_rotation_z_max - @sem_rotation_z_min)
      idx += 1
    end
	
	# freedom 4: case 2
    # lsat variable of snap_to_guide_planes if necessary, with a 50% possibility to resist snapping
    if @model.get_guide_planes.length > 0 && @manifolds <= 1 && !@sem_ignore_guide_planes && x[x.length - 1] < @RATIO_TO_GUIDES
      gplanes = @model.get_guide_planes
      if @manifolds <= 1
        # when not manifolds
        gplane_idx = ([x[x.length - 1].abs / @RATIO_TO_GUIDES, 0.99999999].min * gplanes.length).floor
        # estimated center
        c_x = mx + @offset_x + sx * @len_x_org / 2
        c_y = my + @offset_y + sy * @len_y_org / 2
        c_z = mz + @offset_z + sz * @len_z_org / 2
        # snap to idx
        snap_x, snap_y, snap_z = gplanes.at(gplane_idx).vector_to_point(c_x, c_y, c_z)
        mx -= snap_x
        my -= snap_y
        mz -= snap_z
      else
        # Nothing
      end
    end
    return mx + @offset_x, my + @offset_y, mz + @offset_z, sx, sy, sz, rz
  end
  def bounds(forHorizontal = false)
    if forHorizontal
      # for windows, doors, etc.
      mybox = @instance.bounds()
      bbbox = @base.bounds(false)
      basemin = mybox.min
      basemax = mybox.max
      if (@sem_scale_x_max == @sem_scale_x_min)
        basemin.x = basemax.x
      end
      if (@sem_scale_y_max == @sem_scale_y_min)
        basemin.y = basemax.y
      end
      basemin.x = [basemin.x, bbbox.min.x].min
      basemax.x = [basemax.x, bbbox.max.x].max
      boundingbox = Geom::BoundingBox.new
      boundingbox.add(basemin, basemax)
      return boundingbox
    else
      # for above walls
      basebox = @base.bounds(false)
      basemin = basebox.min
      basemax = basebox.max
      mybox = @instance.bounds()
      if (@sem_scale_x_max == @sem_scale_x_min)
        basemin.x = mybox.center().x
        basemax.x = mybox.center().x
      end
      if (@sem_scale_y_max == @sem_scale_y_min)
        basemin.y = mybox.center().y
        basemax.y = mybox.center().y
      end
      basemin.z = mybox.max().z
      basemax.z = mybox.max().z
      boundingbox = Geom::BoundingBox.new
      boundingbox.add(basemin, basemax)
      return boundingbox
    end 
  end
  def isVerticalRelation(relation)
    return ((relation == "ABOVE" || relation == "BELOW") ? true : false)
  end
  def set_visible(vis)
    @instances.each do |item|
      item.visible = vis
    end
  end
  def update_parent_x()
	if @manifolds > 1
      mymax = (-1E+6)
	  mymin = (1E+6)
	  idx = 0
      while idx < @instances.length do
        mymax = [mymax, @instances.at(idx).bounds.max.x].max
        mymin = [mymin, @instances.at(idx).bounds.min.x].min
        idx +=1
      end
	  @base.ensure_x(mymin, mymax)
	end
  end
  def ensure_x(vmin, vmax)
    if vmin < @gmx || vmax >  @gmx + @gsx * @len_x_org
      mymin = [@gmx, vmin].min
	  mymax = [@gmx + @gsx * @len_x_org, vmax].max
	  @gmx = mymin
	  @gsx = (mymax - mymin) / @len_x_org
	  rollback()
	end
  end
end

class Omg_components
  def initialize
    @array = Array.new
    @num_variables = 0
    @last_component = nil
	@adapter = nil
  end
  def get_last_component()
    return @last_component
  end
  def set_adapter(adp)
    @adapter = adp
  end
  def get_adapter()
    return @adapter
  end
  def add_instance(inst)
    @array.push(inst)
    #if inst.get_num_variables() > @num_variables
      @num_variables = inst.get_num_variables()
	#puts "num_var = #{@num_variables}"
    #end
  end
  def get_ID(x)
    idx = 0
    if @array.length == 1
      return 0
    elsif @array.length > 1
      return (x[@num_variables]*@array.length).floor % @array.length
    end
    return nil
  end
  def record_instance(x, finalize = nil)
    finalize ||= false
    if @array.length > 1
      idx = get_ID(x)
      @array.each do |item|
        # hide else
        item.set_visible(false)
      end
      @array.at(idx).set_visible(true)
      @last_component = @array.at(idx)
    else
      @last_component = @array.at(0)
    end
    if finalize
      if @last_component.get_force_using_guides()
        gplane_idx = ([x[x.length - 1], 0.999999999].min * @last_component.get_model().get_guide_planes().length).floor
        @last_component.get_model().remove_guide_plane(gplane_idx)
        #puts "Removing guide_planes #{gplane_idx}"
      end
      @last_component.record_params()
    end
  end
  def record_best_param()
    @last_component = @array.at(0) if @last_component.nil?
    @last_component.record_best_param()
  end
  def rollback()
    @last_component.rollback()
	show_all()
  end
  def set_best_solution(x)
    @best_solution = x
  end
  def get_force_using_guides()
    return false if @array.length != 1 || @array.at(0).get_manifolds() > 1
    return @array.at(0).get_force_using_guides()
  end
  def transform(x)
    idx = get_ID(x)
    @array.at(idx).transform(x)
  end
  def get_num_variables()
    #if @array.length > 1
    #  return @num_variables + 1
    #end
    return @num_variables
  end
  def reload_num_variables()
    @num_variables = @last_component.get_num_variables()
  end
  def get_guide_planes()
    if @array.length > 1
      # Not applicable for non-manifolds
      return nil
    else
      return @last_component.get_guide_planes()
    end
  end

  def hide_all()
    @array.each do |item|
  	  item.set_visible(false)
    end
  end
  def show_all()
    @array.each do |item|
  	  item.set_visible(true)
    end
  end  
  def update_parent_x()
    @last_component.update_parent_x()
  end
end

class Omg_adapter
  def initialize(components, omg_model)  
    # Instance variables
    @omg_model = omg_model
    @outputfile = @omg_model.get_temp_file()
    @components = components
    @view = Sketchup.active_model.active_view
    @solver = omg_model.get_solver()
    @web_dialog = omg_model.get_web_dialog()
    @iter = 0
    @outputkeys = {
      :filename => @outputfile,
      :width => 640,
      :height => 360,
      :antialias => false,
      :transparent => false
    }
    @best_so_far_adapter = evaluate_fitness()
  end
  def project_to_file()
    @view.write_image @outputkeys
  end
  def evolve(max_iterations)
    @omg_model.set_current_components(@components)
    #@solver.evolve_cmaespp(self, @components.get_num_variables(), max_iterations, @outputfile, "get_measurement", "set_iteration_info")
    @solver.evolve_libcmaes(self, @components.get_num_variables(), max_iterations, @outputfile, "get_measurement", "set_iteration_info", 11)
    #@solver.evolve_galib(self, @components.get_num_variables(), max_iterations, @outputfile, "get_measurement", "set_iteration_info", 0.01, 0.6)
    return @solver.get_best_fitness()
  end
  def evaluate_fitness()
    return @solver.evaluate_pop_representation(@outputfile)
  end
  def get_measurement(x, finalize = nil)
    finalize ||= false
    if (@omg_model.get_guide_planes().length > 1 && @components.get_force_using_guides())
      #puts "Go batching"
      bestf = 1.0
      bestidx = 0
      # batch n times
      (0..(@omg_model.get_guide_planes().length - 1)).each do |idx|
        x[x.length-1] = (idx + 0.5) / @omg_model.get_guide_planes().length
        @components.transform(x)
        #@components.record_instance(x, finalize)
        @view.refresh
        UI.refresh_inspectors()
	    project_to_file()
        f = evaluate_fitness()
        @omg_model.set_iteration_info(f)
        #puts f
        if f < bestf
          bestf = f
          bestidx = idx
        end
      end
      x[x.length-1] = (bestidx + 0.5) / @omg_model.get_guide_planes().length
    end
    @components.transform(x)
    @components.record_instance(x, finalize)
    @view.refresh
    UI.refresh_inspectors()
	project_to_file()
  end
  def set_iteration_info(value)
	@iter += 1
    if (@omg_model.get_guide_planes().length < 2 || !@components.get_force_using_guides())
      @omg_model.set_iteration_info(value)
    end
    if (@best_so_far_adapter > value)
      @best_so_far_adapter = value
      begin
        FileUtils.cp_r @outputfile, @omg_model.get_best_file(), :remove_destination => true
      rescue => detail
        puts "Error running script: " + detail.backtrace.join("\n")
        puts "OMG> Ignored a best-so-far file due to error writing"
      end
    end
    if @omg_model.get_verbose()
      puts  "#{@iter}: #{value}"
    end
  end
  def place_component_from_solution(is_final)
    @best_solution = @solver.ruby_get_best_solution()
    get_measurement(@solver.ruby_get_best_solution(), is_final)
	f = evaluate_fitness()
	FileUtils.cp_r @outputfile, $omg_temp_folder + "best" + f.to_s + ".bmp", :remove_destination => true
    o_time = @solver.get_overall_time().round(4)
    c_time = @solver.get_client_time().round(4)
    e_time = @solver.get_evaluation_time().round(4)
    a_time = @solver.get_algorithm_time().round(4)
    @omg_model.add_time(o_time, c_time, e_time, a_time)
    if @omg_model.get_verbose()
      print "Best fitness: " + (1.0-@solver.get_best_fitness()).round(6).to_s + "\n"
      print "In overall " + o_time.to_s + "s: Client " + c_time.to_s + "s, evaluation " + e_time.to_s + "s, algorithm " + a_time.to_s + "s\n"
    end
  end
  def place_component_from_best_solution()
    get_measurement(@best_solution, true)
	@components.set_best_solution(@best_solution)
    if @omg_model.get_verbose()
      puts "Rollback to best solution..." 
    end
  end
end


def omg_main_proc()
  omg_model = nil
  wd=UI::WebDialog.new( "OMG (Optimization-baed Model Generator) Toolbar", true, "", 
	600, 300, 0, 0, false )
    
  # clear content, reset model
  wd.add_action_callback("omg_reset_model") do |js_wd, omg_target_image|
    Sketchup.active_model.active_entities.clear!()
    $omg_ground = Omg_ground.new(-410, 410, -25, 25, 0)
    omg_model = Omg_model.new
    omg_model.set_geolocation(22.283343, 114.136898)
    omg_model.set_sunrise_time(1980, 6, 21, 0)
	if not omg_target_image.include? ":"
		omg_target_image = File.dirname(__FILE__) + "/omg_modeler/" + omg_target_image
	end
    omg_model.set_ref_measurement(omg_target_image)
    omg_model.set_web_dialog(wd)
    
    # reset cam
    model = Sketchup.active_model
    eye = [0,850,157.5]
    target = [0,0,157.5]
    up = [0,0,1]
    my_camera = Sketchup::Camera.new eye, target, up
    model.active_view.camera = my_camera
    model.active_view.refresh
    UI.refresh_inspectors()
  end
  
  # reset cam
  wd.add_action_callback("omg_reset_cam") do |js_wd|
    # reset view
    model = Sketchup.active_model
    eye = [0,850,157.5]
    target = [0,0,157.5]
    up = [0,0,1]
    my_camera = Sketchup::Camera.new eye, target, up
    view = model.active_view
    view.camera = my_camera
    view.refresh
    UI.refresh_inspectors()
	
  end
  
  # incremental build
  wd.add_action_callback("omg_incremental_apriori") do |js_wd, omg_iters_per_component|
    omg_model.set_iterations_incremental(omg_iters_per_component.to_i)
	omg_model.reset_guide_planes()
    # add a wall on ground (G/F) (2 candidates)
    $omg_ground = Omg_ground.new(-433, 433, -39, 39, 0)
    #omg_model.set_iterations_incremental(omg_iters_per_component.to_i*6)
    gf_facades = omg_model.build_incremental("Components/omg/wall/", ["bars.skp"], "ABOVE", $omg_ground)
    #standing_box = gf_facades.get_last_component().bounds(true)
	#mybox = gf_facades.get_last_component().bounds()
	#UI::messagebox "#{standing_box.min} #{standing_box.max} #{mybox.min} #{mybox.max}"
    # add a wall on 1/F (2 candidates)
    sleep(0.5)
    facades_1f = omg_model.build_incremental("Components/omg/wall/", ["white.skp"], "ABOVE", gf_facades)
	#mybox = facades_1f.get_last_component().bounds()
	#UI::messagebox "#{standing_box.min} #{standing_box.max} #{mybox.min} #{mybox.max}"
  
	
    # add 5 x3section windows to 1/F
    sleep(0.5)
    #omg_model.set_iterations_incremental(omg_iters_per_component.to_i)
    windows_1f = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , facades_1f, 5)
    omg_model.add_guide_planes(windows_1f.get_guide_planes())
	windows_1f.update_parent_x()
	
    #win1f_1 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , facades_1f)
    #win1f_2 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , facades_1f)
    #win1f_3 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , facades_1f)
    #win1f_4 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , facades_1f)
    #win1f_5 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , facades_1f)
  
	
    # add a door on G/F
    sleep(0.5)
    omg_model.set_iterations_incremental(((omg_iters_per_component.to_i-1) / omg_model.get_guide_planes().length).round + 1)
    gf_door = omg_model.build_incremental("Components/omg/door/", ["door1.skp"], "CONTAINS_ON" , gf_facades, nil, true)
	
	
    # add a tree on foreground (2 candidates)
    sleep(0.5)
	omg_model.set_iterations_incremental(omg_iters_per_component.to_i)
    standing_box = gf_facades.get_last_component().bounds(true)
    fore_ground = Omg_ground.new(-350, 350, standing_box.max.y.to_f+20, standing_box.max.y.to_f+220, 0)
    gf_tree = omg_model.build_incremental("Components/omg/tree/", ["palm.skp", "oak.skp"], "ABOVE", fore_ground)
	
	
    # add 3 traditional windows to G/F
    sleep(0.5)
    omg_model.set_iterations_incremental(((omg_iters_per_component.to_i-1) / omg_model.get_guide_planes().length).round + 1)
    windows_Gf1 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , gf_facades, nil, true)
    omg_model.set_iterations_incremental(((omg_iters_per_component.to_i-1) / omg_model.get_guide_planes().length).round + 1)
    windows_Gf2 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , gf_facades, nil, true)
    omg_model.set_iterations_incremental(((omg_iters_per_component.to_i-1) / omg_model.get_guide_planes().length).round + 1)
    windows_Gf3 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , gf_facades, nil, true)
    omg_model.set_iterations_incremental(((omg_iters_per_component.to_i-1) / omg_model.get_guide_planes().length).round + 1)
    windows_Gf4 = omg_model.build_incremental("Components/omg/window/", ["traditional.skp", "french3section.skp"], "CONTAINS_ON" , gf_facades, nil, true)
    
	omg_model.set_refine_mode(true)
	
    #omg_model.redraw_texture()
    omg_model.set_iterations_incremental(omg_iters_per_component.to_i)
    objs = Array.new
	objs << facades_1f
	objs << windows_Gf1
	objs << windows_Gf2
	objs << windows_Gf3
	objs << windows_Gf4
	objs << gf_tree
	objs << gf_door
	omg_model.build_refine(objs)
    omg_model.redraw_texture()
  end
  
  
  # discovery build
  wd.add_action_callback("omg_discovery") do |js_wd, omg_iters_per_component|
    omg_model.set_iterations_incremental(omg_iters_per_component.to_i)
    $omg_ground = Omg_ground.new(-410, 410, -25, 150, 0)
    # discover an obj on ground (G/F) (4 candidates)
    gf_obj = omg_model.build_incremental("Components/omg/", ["tree/oak.skp", "tree/palm.skp", "wall/bars.skp", "wall/white.skp"], "ABOVE", $omg_ground)
    gf_obj1 = omg_model.build_incremental("Components/omg/", ["tree/palm.skp", "tree/oak.skp"], "ABOVE", $omg_ground)
    
    omg_model.redraw_texture()
  end
  
  # export on Google Earth
  wd.add_action_callback("omg_export_kmz") do |js_wd|
    Sketchup.active_model.export $omg_temp_folder + 'sample.kmz'
    system("start " + $omg_temp_folder +  + "sample.kmz")
  end
  
  # load html to GUI
  wd.set_file( File.dirname(__FILE__) + "/omg_modeler/interface.html" )
  wd.show()
  system("explorer file:///" + File.dirname(__FILE__) + "/omg_modeler/monitor.html")
end

UI.menu("Extensions").add_item("OMG: Optimization-baed Model Generator") {
  omg_main_proc()
}
def rotate(axis, angle)
    sel = Sketchup.active_model.selection[0]
    tr = Geom::Transformation.rotation(sel.bounds.center, axis, angle.degrees)
    sel.transform!(tr)
end
