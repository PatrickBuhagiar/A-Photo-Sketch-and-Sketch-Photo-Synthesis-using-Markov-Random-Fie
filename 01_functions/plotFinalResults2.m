function plotFinalResults2(ground, fitted)

  close all
  
  figure(1)
  hold on
  set(gca, 'FontSize', 0.05)
  set(gca, 'FontWeight', 'bold')
  ytick = 0:0.2:1;
  set(gca, 'ytick', ytick);
  limit = 0.08;
  xtick = 0:0.02:limit; 
  set(gca, 'xtick', xtick);
  ylabel('Proportion of Images', 'Interpreter', 'tex', 'fontsize', 15)
  xlabel('Normalized Point to Point RMSE', 'Interpreter', 'tex', 'fontsize', 15)
  title('Multi-PIE frontal expressions (no neutral)','Interpreter', 'tex', 'fontsize', 15)
  grid on
  
  
  %% AOMs-PO
  
  n_test_data = size(fitted, 3);
  shape_RMSE = zeros(n_test_data, 1);
 
  for i = 1:n_test_data
    temp = ground(:,:,i);
    face_size = mean(max(temp) - min(temp) + [1 1]);
    shape_RMSE(i) = mean(sqrt(sum((temp - fitted(:,:,i)).^2,2))) / face_size;
  end
  
  list = 0:0.0005:limit;
  cumulative_shape_RMSE = zeros(size(list));
  count = 1;
  for i = list 
    cumulative_shape_RMSE(count) = length(find(shape_RMSE < i)) / n_test_data;
    count = count + 1;
  end
  
  plot(list, cumulative_shape_RMSE, 'green-^', 'MarkerSize', 5, 'linewidth', 1);
  legend('AOMs - PO', 'Location', 'SouthEast')
  
  
end