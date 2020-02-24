

abstract type abstract_simulation end

include("../examples/ysz_experiments.jl")


function TCtoT(TC)
  return TC+273.15
end

function pO2tosim(pO2)
  return (pO2/100.0 + 1.0e-5)
end



function filename_format_prms(; save_dir="./nouze/", prefix="", prms=Nothing, prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple)

  function consistency_check()
    if (size(prms_names,1) != size(scripted_tuple,1) || 
        size(prms_names,1) != size(prms,1))
      return false
    end
    return true
  end
  
  if !consistency_check()
    println("ERROR: file_format_prms(): shape mismatch (!consistency_check())")
    return throw(Exception)
  end

  scripted_dir = ""
  out_name = prefix
  for i in 1:size(scripted_tuple, 1)
    if scripted_tuple[i] == 1
      scripted_dir = string(scripted_dir,"_",prms_names[i],@sprintf("%.2f",prms[i]))
    end
    out_name = string(out_name, "_", prms_names[i], @sprintf("%.2f",prms[i]))
  end
  
  out_path = string(save_dir,scripted_dir,"/")
  out_name = string(out_name, ".csv")

  (out_path, out_name)
end

function is_between(x, a, b)
    if (b - a) > 0
        if (a <= x) & (x <= b)
            return true
        end
    else
        if (b <= x) & (x <= a)
            return true
        end
    end
    return false
end
