function Movie_Maker_1d(movie_plot,X,no_of_plots,file_name)

title = strcat(file_name,'.avi');

writerObj = VideoWriter(title);
open(writerObj);

plot(X,squeeze(movie_plot(1,:)))
axis([min(X) max(X) 0 1.25]) % Note you can adjust the vertical size of the window by changing the numeric value. 

set(gca,'nextplot','replacechildren');

for mm = 1:no_of_plots
    
    plot(X,squeeze(movie_plot(mm,:)))
    axis([min(X) max(X) 0 1.25]) % Note you can adjust the vertical size of the window by changing the numeric value. 

    frame = getframe;
    writeVideo(writerObj,frame);

end

close(writerObj);
