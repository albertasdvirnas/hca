function [fullTable,barfragq,barfragr] = create_full_table(res_table,bar1,bar2,calc)
    % create full table where are no "loop overs"
    % so that its easier to calculate true positives, false negatives
    
    % also need to take into account when a super-long sequence is mapped,
    % that loops around.., then the match table does not provide a correct
    % result. since 218:218 might mean looping around
    if nargin< 4
        calc = 0;
    end
    
    %todo: check that bar1 and bar2 has same dimensions (rows or cols)
    fullTable = [];
    N = length(bar1);
    M = length(bar2);
    
    barfragq = cell(1,size(res_table, 1));
    barfragr =  cell(1,size(res_table, 1));
    % [res_table, ~] = parse_vtrace(vitResults);
    % problem: both bar1 and bar2 can be circularly shifted..
    for i = 1:size(res_table, 1)
        if res_table(i,1)==0
            res_table(i,1) = N;
        end    
        if res_table(i,2)==0
            res_table(i,2) = N;
        end    
        if res_table(i,3)==0
            res_table(i,3) = M;
        end   
        if res_table(i,4)==0
            res_table(i,4) = M;
        end   
        
        A = res_table(i,1);
        B = res_table(i,2);
        C = res_table(i,3);
        D = res_table(i,4);
        
        if res_table(i,5)==1
            if A > B
                if C > D
                    if N-A == M-C
                        tempTable = ones(2,5);
                        tempTable(1,1:4) = [A N C M];
                        tempTable(2,1:4) = [1 B 1 D ];
                    else
                        tempTable = ones(3,5);
                        if N-A > M-C
                            st = M-C;
                            tempTable(1,1:4) = [A A+st C M];
                            tempTable(2,1:4) = [A+st+1 N 1 N-A-st];
                            tempTable(3,1:4) = [1 B N-A-st+1 D];
                        else
                            st = N-A;
                            tempTable(1,1:4) = [A N C C+st];
                            tempTable(2,1:4) = [1 M-C-st C+st+1 M];
                            tempTable(3,1:4) = [ M-C-st+1 B 1 D];
                        end
                    end
                else
                    tempTable = ones(2,5);
                    st = N-A;
                    tempTable(1,1:4) = [A N C C+st];
                    tempTable(2,1:4) = [1 B C+st+1 D];
                end
            else
                if C > D
                    tempTable = ones(2,5);
                    st = M-C;
                    tempTable(1,1:4) = [A A+st C M];
                    tempTable(2,1:4) = [A+st+1 B 1 D];
                else
                    tempTable = ones(1,5);
                    tempTable(1:4) = res_table(i,1:4);
                end
            end
        else
            if A > B
                if C < D
                    if N-A+1 == C
                        tempTable = ones(2,5);
                        tempTable(1,1:4) = [A N C 1];
                        tempTable(2,1:4) = [1 B M D];
                    else
                        tempTable = ones(3,5);
                        if N-A+1 > C
                            st = C;
                            tempTable(1,1:4) = [A A+st-1 C 1];
                            tempTable(2,1:4) = [A+st N M M-(N-A-st+1)+1]; % maybe +1?
                            tempTable(3,1:4) = [1 B M-(N-A-st+1) D];
                        else
                            st = N-A;
                            tempTable(1,1:4) = [A N C C-st];
                            tempTable(2,1:4) = [1 C-st-1 C-st-1 1];
                            tempTable(3,1:4) = [C-st B M D];
                        end
                    end
                else
                    tempTable = ones(2,5);
                    st = N-A;
                    tempTable(1,1:4) = [A N C C-st];
                    tempTable(2,1:4) = [1 B C-st-1 D];
                end
            else
                if C < D
                    tempTable = ones(2,5);
                    st = C;
                    tempTable(1,1:4) = [A A+st-1 C 1];
                    tempTable(2,1:4) = [A+st B M D];
                else
                    tempTable = ones(1,5);
                    tempTable(1:4) = res_table(i,1:4);
                end
            end            
            
            
            tempTable(:,5) = 2;

        end
            
%         for j=1:size(tempTable,1)
%             if res_table(i,5)==2  % only second sequence can be inverted
%                 tempTable(:,5) = 2;
%             end
% 
%         end
%             

        
        fullTable = [ fullTable;  tempTable];
        if calc == 1
            barfragq{i}= [];
            barfragr{i} = [];
            for j=1:size(tempTable,1)
                % here check if we're not taking too few..
                %
                %
                 barfragq{i} = [ barfragq{i} bar1(tempTable(j,1):tempTable(j,2))];
                 if tempTable(j,5) == 1
                      barfragr{i} =[ barfragr{i}; bar2(tempTable(j,3):tempTable(j,4))'];
                 else
                     barfragr{i} =[ barfragr{i}; bar2(tempTable(j,3):-1:tempTable(j,4))'];
                 end
                 lenDiff  = length(barfragr{i}) -length(barfragq{i});
                  if abs(lenDiff) > 0
                     if lenDiff > 0
                         if tempTable(j,2)+lenDiff > length(bar1)
                            barfragq{i} = [  barfragq{i} bar1(tempTable(j,2)+1:end)  bar1(1:lenDiff-length(bar1)+tempTable(j,2))];
                         else
                             barfragq{i} = [  barfragq{i} bar1(tempTable(j,2)+1:tempTable(j,2)+lenDiff)];
                         end
                     else
                        if tempTable(j,4)+abs(lenDiff) > length(bar2)
                            barfragr{i} = [  barfragr{i} bar2(tempTable(j,4)+1:end)'  bar2(1:abs(lenDiff)-length(bar2)+tempTable(j,4))'];
                         else
                             barfragq{i} = [  barfragr{i} bar2(tempTable(j,4)+1:tempTable(j,4)+abs(lenDiff))'];
                         end 
                     end
                  end  
            end
        end

%                    fullTable = [ fullTable;  tempTable];
%         if plot == 1
%         end
    end
    
end

