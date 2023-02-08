function bits_hat = ldpc_bp_function(L, H, max_iter) %#codegen
nb_c_nodes = size(H,1); % the number of check nodes (= number of rows)
nb_v_nodes = size(H,2); % the number of v_nodes (= number of columns)

for idx_iter_llr = 1:1:max_iter
                % 1) Initial LLRs from VNs to CNs
                if(idx_iter_llr == 1)
                    M = H .* L';
                end
                
                % 2) LLRs from CNs to VNs
                E = zeros(nb_c_nodes, nb_v_nodes);
                for idx_row = 1:1:nb_c_nodes
                    m_tmp = M(idx_row, :);
                    [~, col_tmp] = find(H(idx_row, :));
                    m_tmp = m_tmp(:,col_tmp);
                    m_tmp(m_tmp < -4) = -4;
                    m_tmp(m_tmp > 4) = 4;
                    matrix = eye(length(m_tmp)) < 1;
                    matrix_tmp = matrix .* m_tmp;
                    m_tmp_1 = tanh(matrix_tmp/2);
%                     m_tmp_2 = atanh(prod(m_tmp_1 + eye(length(m_tmp)), 2));
                    x = prod(m_tmp_1 + eye(length(m_tmp)), 2);
                    m_tmp_2 = log((1 + x)./(1 - x));
                    E(idx_row, col_tmp) = 2 * m_tmp_2;
%                     for idx_col = 1:1:length(col_tmp)
%                         tmp = m_tmp;
%                         tmp(:,idx_col) = 1;
%                         res = atanh(prod(tmp));
%                         res(res < -19.07) = -19.07;
%                         res(res > 19.07) = 19.07;
%                         E(idx_row, col_tmp(idx_col)) = 2*res;
%                     end
                end

                % 3) LLRs from VNs to CNs
                M = zeros(nb_c_nodes, nb_v_nodes);
                for idx_col = 1:1:nb_v_nodes
                    e_tmp = E(:, idx_col);
                    [row_tmp, ~] = find(H(:, idx_col));
                    e_tmp = e_tmp(row_tmp, :);
                    if(length(e_tmp) ~= 1)
                        matrix = eye(length(e_tmp)) < 1;
                        matrix_tmp = matrix .* e_tmp;
                        e_tmp_2 = sum(matrix_tmp, 1);
%                         e_tmp_2(e_tmp_2 < -19.07) = -19.07;
%                         e_tmp_2(e_tmp_2 > 19.07) = 19.07;
                        M(row_tmp, idx_col) = e_tmp_2;
                    else 
                        M(row_tmp, idx_col) = e_tmp;
                    end
                end

                % Update total LLRs
                E_tmp = sum(E, 1);
                L = L + E_tmp';
                
                % Check if a correct codeword is found 
                bits_hat = L < 0;
                if(any(mod(H * bits_hat, 2)) == 0)
                    break;
                end
end   
end
          