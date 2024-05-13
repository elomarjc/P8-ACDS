function out = skew_matrix(x)
    assert(size(x,1) == 3, 'skew_matrix: Input argument is of wrong dimension (not 1x3)');
    
    out = [0    -x(3)  x(2);
           x(3)  0    -x(1);
          -x(2)  x(1)  0];
end

