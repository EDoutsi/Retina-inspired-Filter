function F = Matrix_reshape(F_Matrix,i,Fx,Fy)
     F = F_Matrix(i,:);
     F = reshape(F,Fx,Fy);
end