% Sorting Data
% Author: Raphael Gilioli Heineck, Chemical Eng.

function sortData

clear;
clc;
close all;

par = xlsread('adsorptionData.xlsx','dadosIsoterma','X6:Z14');
cexqe = adsorptionCapacity(par);

ce = cexqe(:,1);
qe = cexqe(:,2);
n = length(ce);

ceZ = zeros(n*2,1);
qeZ = zeros(n*2,1);

data = zeros(n*2,1);

for n = 1:n
    if n == 1
        ceZ(n) = ce(n);
        qeZ(n+1) = qe(n);
    end
end



% for n = 0:2:n
%     j = rem(n,2);
%     if j == 1
%         ceZ(n) = ce(n-1);
%     end
% end
% 
% for n = 1:n
%     j = rem(n,2);
%     if j == 0
%         qeZ(n) = qe(n-1);
%     end
% end



% for n = 1:n
%     j = rem(n,2);
%     if j == 0
%         qeZ(n) = qe(n-1);
%     else
%         if j == 1
%             ceZ(n) = ce(n-1);
%         end
%     end
% end
% 
% for n = 1:n
%     j = rem(n,2);
%     if j == 0
%         data(n) = qeZ(n);
%     else
%         if j == 1
%             data(n) = ceZ(n);
%         end
%     end
% end
% 
% disp(data);

% for n = 1:2:n
%     data(n) = ce(n);
%     data(n+1) = qe(n);
% end

% for n = 1:n
%     data(n) = [ce(n) qe(n)];
% end

% for n = 1:3:n
%     data(n) = cexqe(n,m-1);
%     data(n+1) = cexqe(n+1,m);
% end


end