function sum = calc_ISI(P,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the ISI given a P matrix, and a length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sum_1 = 0;
        sum_2 = 0;
        sum_1a = 0;
        sum_2a = 0;
        sum = 0;
        for i = 1:L
            for j = 1:L
               sum_1 = sum_1 + ((abs(P(i,j))/max(abs(P(i,:)))));
            end
            sum_1a = sum_1a + sum_1 - 1;
            sum_1 = 0;
        end

        for j = 1:L
            for i = 1:L
               sum_2 = sum_2 + ((abs(P(i,j))/max(abs(P(:,j)))));
            end
            sum_2a = sum_2a + sum_2 - 1;
            sum_2 = 0;
        end

        sum = (1/(2*L*(L - 1)))*(sum_1a + sum_2a);
    end
