function bits_hat = ldpc_bp_function(L, H, max_iter) %#codegen
nb_c_nodes = size(H,1); % the number of check nodes (= number of rows)
nb_v_nodes = size(H,2); % the number of v_nodes (= number of columns)

% 1) Initial LLRs from VNs to CNs
M = H .* L';
L_int = L;

for idx_iter_llr = 1:1:max_iter             
                % 2) LLRs from CNs to VNs
                E = zeros(nb_c_nodes, nb_v_nodes);
                for idx_row = 1:1:nb_c_nodes
                    h_tmp = logical(H(idx_row,:));
                    m_tmp = M(idx_row, :);
                    m_tmp = m_tmp(:, h_tmp);
                    matrix = eye(length(m_tmp)) < 1;
                    matrix_tmp = matrix .* tanh(m_tmp/2);
%                     m_tmp_1 = matrix_tmp;
%                     m_tmp_2 = 2*atanh(prod(m_tmp_1 + eye(length(m_tmp)), 2));
                    x = prod(matrix_tmp + eye(length(m_tmp)), 2);
%                     m_tmp_2 = log((1 + x)./(1 - x));
                    E(idx_row, h_tmp) = log((1 + x)./(1 - x));
%                     E(idx_row, h_tmp) = 2*atanh(x);
                end

                % 3) LLRs from VNs to CNs
                M = zeros(nb_c_nodes, nb_v_nodes);
                for idx_col = 1:1:nb_v_nodes
                    h_tmp = logical(H(:, idx_col));
                    e_tmp = E(h_tmp, idx_col);
%                     e_tmp = e_tmp(h_tmp, :);
                    if(length(e_tmp) ~= 1)
                        matrix = eye(length(e_tmp)) < 1;
                        matrix_tmp = matrix .* e_tmp;
                        e_tmp_2 = sum(matrix_tmp, 1);
                        M(h_tmp, idx_col) = e_tmp_2 + L_int(idx_col) ;
                    else 
                        M(h_tmp, idx_col) = e_tmp + L_int(idx_col);
                    end
                end

                % Update total LLRs
                E_tmp = sum(E, 1);
                L = L_int + E_tmp';
                
                % Check if a correct codeword is found 
                bits_hat = L < 0;
                if(any(mod(H * bits_hat, 2)) == 0)
                    break;
                end
end   
end
          
