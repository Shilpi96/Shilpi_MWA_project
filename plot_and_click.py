def onclick(event):
    global ix
    ix= event.xdata
    global iy
    iy=event.ydata
    global x_coords
    x_coords.append(ix)
    global y_coords
    y_coords.append(iy)
    #print(len(x_coords))
    # Disconnect after 2 clicks
    if len(x_coords) == 7:             ####### wanted the x coordinates for 7 points
        fig.canvas.mpl_disconnect(cid)   ######## once you click 7 times the plot will close, 
                                           ##### then x_coords should have all the values you clicked on.
        plt.close()
    return
....
PLOT AS NORMAL HERE
.....
x_coords=[]
y_coords=[]
cid = fig.canvas.mpl_connect('button_press_event', onclick) #2.5 arcsec per pixel
plt.show()

######  to save x_coords to have for the future
# my_saved_values= np.load('x_cood_save.npy')