% GNU Octave code
clear all

function [q_count, q_values, I_values] = read_data ()
  filename = "output_monot_log.txt";
  fid = fopen(filename, "r");
  fscanf(fid, "%s", 1);  # read "q"
  [q_values, q_count, errmsg] = fscanf(fid, "%f", Inf);  # read values of q
  fscanf(fid, "%s", 1);  # read "I(q)"

  I_values = [];
  while (1)
    [val, count, errmsg] = fscanf(fid, "%s", 3);  # read "m = value"
    if (length(val) == 0)
      break;
    endif
    [val, count, errmsg] = fscanf(fid, "%f", q_count);  # read values of I(q)
    I_values = [I_values; val'];
  endwhile
  fclose(fid);
endfunction


[q_count, q_values, I_values] = read_data();

while (1)
  m = input("m = ");  # time step to plot
  if (m == 0)
    break;
  endif
  plot(q_values, I_values(m,:));
endwhile
