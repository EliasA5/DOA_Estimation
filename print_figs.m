%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: helper function to print all open figures, save them in a folder
% called pics in the same directory as the file, and does not override
% existing photos as long as they are named in the corrent order (1,2,3...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h =  findobj('type','figure');
fig_count = length(h);
res = dir("./pics/*.png");
for i=(1:fig_count)+length(res)
    name = num2str(i);
    saveas(figure(i-length(res)),append('./pics/',name),'png')
end