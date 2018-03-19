function [output , state] = BCH(data )

reg = [0,0,0,0];
state = 0 ;

for i=1 :15
    savebit1 = xor(reg(4),data(i));
    savebit2 = xor(reg(4),reg(1 ));
    reg(3:4) = reg(2:3);
    reg(1)   = savebit1;
    reg(2)   = savebit2;
end
switch num2str(fliplr(reg))
    case '0  0  0  0'
        err_cor = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
        state = 1 ;
    case '0  0  0  1'
        err_cor = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
    case '0  0  1  0'
        err_cor = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0];
    case '0  0  1  1'
        err_cor = [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0];
    case '0  1  0  0'
        err_cor = [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];
    case '0  1  0  1'
        err_cor = [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0];
    case '0  1  1  0'
        err_cor = [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0];
    case '0  1  1  1'
        err_cor = [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0];
    case '1  0  0  0'
        err_cor = [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0];
    case '1  0  0  1'
        err_cor = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    case '1  0  1  0'
        err_cor = [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0];
    case '1  0  1  1'
        err_cor = [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0];
    case '1  1  0  0'
        err_cor = [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0];
    case '1  1  0  1'
        err_cor = [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0];
    case '1  1  1  0'
        err_cor = [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
    case '1  1  1  1'
        err_cor = [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
end
if (~state)
    fprintf('error');
end
output = xor(data , err_cor);
output = output+0;