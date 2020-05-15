x = 1:10;% 10 stimulus
y = 1:8;% 8 subparts
z = randi(5,8,10)-1;% random value from 0 to 4
imagesc(x,y,z)
colormap(jet(5))