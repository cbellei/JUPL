function is_print_time(tm, dt_print)
  if mod(tm, dt_print) < dtm
      is_print = true
      println("--writing files, t = ", tm, "; maxTime = ", maxTime, "; dtm = ", dtm)
  else
      is_print = false
  end
  return is_print
end


function is_quiet_time!(time1, time2, tm, tm_quiet, geom, nz, nq)
    if (time1 && tm < tm_quiet)
        if geom=="slab"
            nq = nquiet
        elseif geom=="spherical"
            nq = nz - nquiet
        end
        println("nz = ", nz, " nquiet = ", nquiet, " tm_quiet=", tm_quiet)
        time1 = false
    elseif (time2 && tm >= tm_quiet)
        println("end of quiet time")
        nq = nz
        # dtm = dtm / dt_multi
        # lm = dtm / dr
        time2 = false
    end
    return nq, time1, time2
end
