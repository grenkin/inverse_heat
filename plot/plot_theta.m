% GNU Octave code
clear all

function [x_count, x_values, theta_values] = read_data ()
  filename = "output_theta.txt";
  fid = fopen(filename, "r");
  fscanf(fid, "%s", 1);  # read "x"
  [x_values, x_count, errmsg] = fscanf(fid, "%f", Inf);  # read values of x
  fscanf(fid, "%s", 2);  # read "theta(x, t_m)"

  theta_values = [];
  while (1)
    [val, count, errmsg] = fscanf(fid, "%s", 3);  # read "m = value"
    if (length(val) == 0)
      break;
    endif
    [val, count, errmsg] = fscanf(fid, "%f", x_count);  # read values of theta(x, t_m)
    theta_values = [theta_values; val'];
  endwhile
  fclose(fid);
endfunction


[x_count, x_values, theta_values] = read_data();

while (1)
  m = input("m = ");  # time step to plot
  if (m == -1)
    break;
  endif
  plot(x_values, theta_values(m + 1, :));
endwhile
