function draw_boundaries(Mdl,type)
% DRAW_BOUNDARIES 
% This function on draws class boundaries based on the model input
% Written by Evan Campbell for ECE4553, 2019

xlims = xlim;
ylims = ylim;

hold on

switch(type)
    case 'LDA'
        num_class = length(Mdl.ClassNames);

        % cycle through the class combinations (in my case only two
        % classes, so only one dividing line)
        for c1 = 1:num_class
            for c2 = 1:num_class
                if c1 == c2
                    continue
                elseif c1 < c2
                    continue
                end
                
                % return parameters for the equation describing the linear
                % division
                K = Mdl.Coeffs(c1,c2).Const;
                L = Mdl.Coeffs(c1,c2).Linear;
                
%                 fprintf('\nf_LDA = @(x1,x2) K + L(1)*x1 + L(2)*x2')
%                 fprintf('\nconstant K: %f', K)
%                 fprintf('\ncoefficient L(1): %f', L(1))
%                 fprintf('\ncoefficient L(2): %f\n', L(2))

                % the equation of the seperating function
                % x1 = x
                % x2 = y
                f_LDA = @(x1,x2) K + L(1)*x1 + L(2)*x2;

                % when f_LDA = 0 the equation can be rearranged to get 
                % y = mx + b. Find m and b when that is done
                b = -K/L(2);
                m = -L(1)/L(2);
                line_eq = sprintf('\ny = %.4fx + %.4f\n', m, b);
                fprintf(line_eq)

                % graph the line function when f_LDA is set to 0
                h3 = fimplicit(f_LDA,[xlims ylims]);
                h3.Color = 'k';

%                 % add equation to the plot
%                 text(xlims(2)-30, ylims(1)+0.5, line_eq)
%                 h3.LineWidth = 2;
    
            end
        end
    case 'QDA'
        num_class = length(Mdl.ClassNames);
        for c1 = 1:num_class
            for c2 = 1:num_class
                if c1 == c2
                    continue
                elseif c1 < c2
                    continue
                end
                K = Mdl.Coeffs(c1,c2).Const;
                L = Mdl.Coeffs(c1,c2).Linear;
                Q = Mdl.Coeffs(c1,c2).Quadratic;
                
                f_QDA = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
                    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
                h3 = fimplicit(f_QDA,[xlims ylims]);
                h3.Color = 'c';
                h3.LineWidth = 2;
            end
        end
    case 'kNN'
        
        % Train KNN
        
        grid_x = xlims(1):0.1:xlims(2);
        grid_y = ylims(1):0.1:ylims(2);
        [xx1,xx2] = meshgrid(grid_x, grid_y);
        XGrid = [xx1(:) xx2(:)];
        predicted_knn = predict(Mdl,XGrid);
        gscatter(xx1(:), xx2(:), predicted_knn,'rbk');
        
end


end

