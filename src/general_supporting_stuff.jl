
mutable struct prms_struct
  names
  values
end

function consistency_check(prms::prms_struct)
  if prms.names == Nothing && prms.values == Nothing
    return true
  end
  
  if size(prms.names,1) != size(prms.values,1)
    println("ERROR: consistency_check: shape mismatch size(prms.names,1) != size(prms.values,1) >>>> $(size(prms.names,1)) != $(size(prms.values,1)) ")
    return throw(Exception)
  end
  
  for i in 1:size(prms.values,1)
    if typeof(prms.values[i]) != Colon
      if size(prms.values[i],1) < 1
        println("ERROR: consistency_check: empty_field prms.values[$i]")
        return throw(Exception)
      end
    end
  end
  
  return true
end

function for_each_prms_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else
      for parameter in prms_lists[active_idx]
        recursive_call(push!(deepcopy(output_set),parameter), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end

function for_each_indicies_in_prms_lists(prms_lists, perform_generic)
  function recursive_call(output_set, active_idx)
    if active_idx > size(prms_lists,1)
      perform_generic(output_set)
    else
      for i in 1:size(prms_lists[active_idx],1)
        recursive_call(push!(deepcopy(output_set),i), active_idx + 1)
      end
    end
  end
  recursive_call([],1)
  return
end

function get_prms_from_indicies(prms_lists, tuple)
  output_array = []
  for (i,list) in enumerate(prms_lists)
    append!(output_array, list[tuple[i]])
  end
  return output_array
end

function check_equal_size(list_1, list_2)
  if length(list_1) == length(list_2)
    return true
  else
    println("ERROR: check_equal_size: shape mismatch $(list_1) and $(list_2)")
    return throw(Exception)
  end
end

function check_x_in(x, low, upp)
  for (i, item) in enumerate(x)
    if item < low[i] || upp[i] < item
      return false
    end
  end
  return true
end


