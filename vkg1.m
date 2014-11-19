%program for 3d plotting.
[X,Y] = meshgrid([1:.5:10],[1:.5:10]);
Z = X .^2 - Y.^2;
surf(X,Y,Z)
figure,mesh(X,Y,Z);

c= cell(2,3)
s = c{1,2}
c{1,1} = [1 2 3;4 5 6]
c{1,2} = 'kiran goud'
c(1,3) = {'vk'}
c(1,4) = {'vicky'}
c{2,1} = 'hello';
c{1,2}(1,1)
c{1,2}(1,:)
