function plotFinalResults4(ground, fitted)

  close all
  
  figure(1)
  hold on
  set(gca, 'FontSize', 0.05)
  set(gca, 'FontWeight', 'bold')
  ytick = 0:0.2:1;
  set(gca, 'ytick', ytick);
  limit = 10;
  xtick = 0:1:limit; 
  set(gca, 'xtick', xtick);
  ylabel('Proportion of Images', 'Interpreter', 'tex', 'fontsize', 15)
  xlabel('Shape RMSE', 'Interpreter', 'tex', 'fontsize', 15)
  title('XM2VTS','Interpreter', 'tex', 'fontsize', 15)
  grid on
  
  
  %% AOMs-PO
  
  n_test_data = size(fitted, 3);
  shape_RMSE = zeros(n_test_data, 1);
 
  for i = 1:n_test_data
    temp1 = ground(:,:,i);
    temp2 = fitted(:,:,i);
    shape_RMSE(i) = sqrt(mean((temp1(:) - temp2(:)).^2));
  end
  
  list = 0:0.1:limit;
  cumulative_shape_RMSE = zeros(size(list));
  count = 1;
  for i = list 
    cumulative_shape_RMSE(count) = length(find(shape_RMSE < i)) / n_test_data;
    count = count + 1;
  end
  
  plot(list, cumulative_shape_RMSE, 'green-^', 'MarkerSize', 5, 'linewidth', 1);
  legend('AOMs - PO', 'Location', 'SouthEast')
  
  
end