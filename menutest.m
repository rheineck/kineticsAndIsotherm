function menutest

p = input('type a number between 1 - 8: ');

while p == 1
    disp('Number One');
    p = input('another number? 1-8 ');
end
while  p == 2
    disp('Number Two');
    p = input('another number? 1-8 ');
end
if p ~= 1 || p~= 2
    disp('end');
end

end