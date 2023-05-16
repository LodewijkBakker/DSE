# Python program to create
# a file explorer in Tkinter

# import all components
# from the tkinter library
from tkinter import *

# import filedialog module
from tkinter import filedialog

import numpy as np
from scipy.spatial.transform import Rotation as R

import os

from project_area_initialiser import torque_calculator
from rot_angles_generators import sphere_fib_user_angle, sphere_fib_octant, pole_accurate_angles


def update_rot_angles(ng, user_input_angles):
    global rot_angles

    octant_sel = [octant_1.get(), octant_2.get(), octant_3.get(),
                  octant_4.get(), octant_5.get(), octant_6.get(),
                  octant_7.get(), octant_8.get()]
    n_octant_sel = np.count_nonzero(octant_sel)  # checks amount of octant selected,
    # n_octant_sel should not be 0 because of errors therefore no octants
    # selected will be all selected as well
    if n_octant_sel == 0:
        octant_sel = [True, True, True, True, True, True, True, True]
        n_octant_sel = 8
    if ng == 0:
        ng = 100  # to prevent errors

    if angle_choice.get() == 2:
        if n_octant_sel != 8:
            ng = int(ng * 4 / n_octant_sel)
            # Update ng for the amount, so that approximately it is corect
            # This leaves allot to be desired
        rot_angles = pole_accurate_angles(ng, octant_sel)
    else:
        ng = int(ng * 8 / n_octant_sel)  # Update ng for the amount, set closest int
        rot_angles = sphere_fib_octant(ng, octant_sel)
        rot_angles = sphere_fib_user_angle(ng)

    draw_angle_point_rep_canvas(rot_angles, user_input_angles)


def browseFiles():
    global filename
    filename = filedialog.askopenfilename(initialdir=os.getcwd(),
                                          title="Select a File",
                                          filetypes=(("Stl files",
                                                      "*.stl*"),
                                                     ("all files",
                                                      "*.*")))
    # Change label contents
    label_file_explorer.configure(text="File Opened: " + filename)


def browse_export_folder():
    global export_folder
    export_folder = filedialog.askdirectory(initialdir=os.getcwd(),
                                          title="Select a folder")
    # Change label contents
    label_folder_export.configure(text="Folder Opened: " + export_folder)


def export_angles(res):
    export_file_path_n_name = export_folder +"/" + export_file_name.get() + ".txt"

    angle_export_file = open(export_file_path_n_name, "w")
    angle_export_file.write("Drag surface m^2 || Area_Arm m^3 || Euler rotx || Euler roty || Euler rotz \n")
    max_drag_surface = max(i[0] for i in res)
    angle_export_file.write("Max Area:  {} \n".format(max_drag_surface))
    max_area_arm = max(i[2] for i in res)  # could maybe be made faster by combining
    # the max above
    angle_export_file.write("Max Area_Arm:  {} \n".format(max_area_arm))  # combo of area and arm max moment
    for i in res:
        single_line_write = str(i[0]) + ", " + str(i[2]) + ", " + str(i[1][0]) + "," + str(i[1][1]) + "," + str(i[1][2]) + "\n"
        # could be done cleaner with fstrings but eh for now it's okay
        angle_export_file.write(single_line_write)


def draw_angle_point_rep_canvas(rot_angles, user_input_angles):
    point_canvas.delete("all")
    point_scale = 85  # how spaced the points are visually
    start_translate_x = 100  # determines where the center x coord is of the points
    start_translate_y = 100  # determines where the center y coord is of the points
    point_size = 1  # how many pixels the diameter is of a drawn point
    text_translate_y = 10  # how much x y and z should be below their respective axis

    rotation_user_input = R.from_euler('zyx', user_input_angles, degrees=True)
    rotation_on_x_vector = rotation_user_input * R.from_euler('zyx', rot_angles, degrees=True)
    projection_points = rotation_on_x_vector.apply([1, 0, 0]) * point_scale
    # applies all rotation that are input by user and by generation

    # Draws 3d lines on object as graph
    x_axis_start = rotation_user_input.apply([-1, 0, 0]) * point_scale
    y_axis_start = rotation_user_input.apply([0, -1, 0]) * point_scale
    z_axis_start = rotation_user_input.apply([0, 0, -1]) * point_scale

    x_axis_end = rotation_user_input.apply([1, 0, 0]) * point_scale
    y_axis_end = rotation_user_input.apply([0, 1, 0]) * point_scale
    z_axis_end = rotation_user_input.apply([0, 0, 1]) * point_scale

    point_canvas.create_line(x_axis_start[0] + start_translate_x,
                             x_axis_start[1] + start_translate_y,
                             x_axis_end[0] + start_translate_x,
                             x_axis_end[1] + start_translate_y,
                             fill='blue')
    point_canvas.create_line(y_axis_start[0] + start_translate_x,
                             y_axis_start[1] + start_translate_y,
                             y_axis_end[0] + start_translate_x,
                             y_axis_end[1] + start_translate_y,
                             fill='red')
    point_canvas.create_line(z_axis_start[0] + start_translate_x,
                             z_axis_start[1] + start_translate_y,
                             z_axis_end[0] + start_translate_y,
                             z_axis_end[1] + start_translate_x,
                             fill='green')

    # dot to show positive side of axis
    point_canvas.create_oval(x_axis_end[0] + start_translate_x - 3,
                             x_axis_end[1] + start_translate_y - 3,
                             x_axis_end[0] + start_translate_y + 3,
                             x_axis_end[1] + start_translate_x + 3,
                             fill='blue')
    point_canvas.create_oval(y_axis_end[0] + start_translate_x - 3,
                             y_axis_end[1] + start_translate_y - 3,
                             y_axis_end[0] + start_translate_y + 3,
                             y_axis_end[1] + start_translate_x + 3,
                             fill='red')
    point_canvas.create_oval(z_axis_end[0] + start_translate_x - 3,
                             z_axis_end[1] + start_translate_y - 3,
                             z_axis_end[0] + start_translate_y + 3,
                             z_axis_end[1] + start_translate_x + 3,
                             fill='green')

    point_canvas.create_text(x_axis_end[0] + start_translate_x,
                             x_axis_end[1] + start_translate_y + text_translate_y, text="+x")

    point_canvas.create_text(y_axis_end[0] + start_translate_x,
                             y_axis_end[1] + start_translate_y + text_translate_y, text="+y")

    point_canvas.create_text(z_axis_end[0] + start_translate_x,
                             z_axis_end[1] + start_translate_y + text_translate_y, text="+z")

    for draw_point in projection_points:
        point_canvas.create_oval(draw_point[0] + start_translate_x,
                                 draw_point[1] + start_translate_y,
                                 draw_point[0] + start_translate_x + point_size,
                                 draw_point[1] + start_translate_y + point_size,
                                 fill="blue",
                                 tag="circ")



if __name__ == "__main__":
    # Create the root window
    window = Tk()

    # Set window title
    window.title('File Explorer')

    # Set window size
    window.geometry("700x500")

    # Set window background color
    window.config(background="white")


    # Create a File Explorer label
    label_file_explorer = Label(window,
                                text="File Explorer using Tkinter",
                                fg="blue", width="50")

    button_explore = Button(window,
                            text="Browse Files",
                            command=browseFiles, width="50")

    # Create Future Folder explorer Label, export buttons and export file name entry
    label_folder_export = Label(window,
                                text="Here path to opened folder can be seen",
                                fg="blue", width="50")
    export_button = Button(window,
                          text="Select place for export", width="50", command=browse_export_folder)

    export_to_folder_button = Button(window,
                          text="Press to export", width="50", command=export_angles)
    export_file_name = StringVar()
    export_file_name_entry = Entry(textvariable=export_file_name, width="20")
    export_file_name_entry_label = Label(window, text="Export file name", width="30")

    button_start = Button(window, text="Calculate Torques",
                          command=lambda: torque_calculator(filename, rot_angles), width="50")

    # Grid method is chosen for placing
    # the widgets at respective positions
    # in a table like structure by
    # specifying rows and columns
    label_file_explorer.grid(sticky="W", columnspan=2, column=0, row=2)

    button_explore.grid(sticky="W", columnspan=2, column=0, row=1)
    button_start.grid(sticky="W", columnspan=2, column=0, row=7)

    angle_choice = IntVar()
    R1 = Radiobutton(window, text="Fibonacci Angles (Accurate Spread)", variable=angle_choice,
                     value=1,
                     width="47", anchor="w")
    R1.grid(sticky="W", columnspan=2, column=0, row=3)

    R2 = Radiobutton(window, text="Polar Circular symmetric Angles", variable=angle_choice,
                     value=2,
                     width="47",
                     anchor="w")
    R2.grid(sticky="W", columnspan=2, column=0, row=4)

    EntryText = Label(window, text="How many points of sim detail approx?", width="30")
    EntryText.grid(sticky="W", columnspan=2, column=0, row=6)

    point_canvas_text = Label(window, text="Points shown", width="30")
    point_canvas_text.grid(sticky="W", column=2, row=1, columnspan="2")

    # Entry widget for amount of sphere points
    n_sphere_points = IntVar()
    n_sphere_points.trace("w", lambda name, index, mode,
                                      var=n_sphere_points: update_rot_angles(n_sphere_points.get(),
                                                                             [user_input_theta.get(),
                                                                              user_input_phi.get(),
                                                                              user_input_psi.get()]))
    n_sphere_points_entry = Entry(width="20", textvariable=n_sphere_points)
    n_sphere_points_entry.grid(sticky="E", columnspan=1, column=1, row=6)

    # # Vector in, both tet and other
    # EntryVecText = Label(window, text="Please enter center of gravity vector from bottom left point", width="50")
    # EntryVecText.grid(sticky="W", columnspan=2, column=0, row=8)
    #
    # grav_center_x = IntVar()
    # grav_center_y = IntVar()
    # grav_center_z = IntVar()
    #
    # grav_center_x_entry = Entry(width="10", textvariable=grav_center_x)
    # grav_center_x_entry.grid(sticky="E", columnspan=1, column=0, row=9)
    #
    # grav_center_y_entry = Entry(width="10", textvariable=grav_center_y)
    # grav_center_y_entry.grid(sticky="E", columnspan=1, column=1, row=9)
    #
    # grav_center_z_entry = Entry(width="10", textvariable=grav_center_z)
    # grav_center_z_entry.grid(sticky="E", columnspan=1, column=2, row=9)

    # export states:
    export_button.grid(sticky="W", columnspan=2, column=0, row=11)
    export_to_folder_button.grid(sticky="W", columnspan=2, column=0, row=14)
    label_folder_export.grid(sticky="W", columnspan=2, column=0, row=12)
    export_file_name_entry_label.grid(sticky="W", column=0, columnspan="2", row=13)
    export_file_name_entry.grid(sticky="E", column=1, row=13)

    # Point canvas
    point_canvas = Canvas(window, width="200", height="200", borderwidth=1, relief="solid")
    point_canvas.grid(column=2, row=2, rowspan="7", columnspan="2")

    slide_horizontal_1_text = Label(window, text="Angle around z (yaw)", width="10")
    slide_horizontal_1_text.grid(column=2, row=9)
    user_input_theta = IntVar()
    user_input_theta.trace("w", lambda name, index, mode, var=user_input_theta: draw_angle_point_rep_canvas(rot_angles, [user_input_theta.get(), user_input_phi.get(), user_input_psi.get()]))
    slide_horizontal_2 = Scale(window, from_=0, to=360, orient=HORIZONTAL,
                               variable=user_input_theta)  # this is not theta
    slide_horizontal_2.grid(column=3, row=9)

    slide_horizontal_2_text = Label(window, text="Angle around y (pitch)", width="10")
    slide_horizontal_2_text.grid(column=2, row=10)
    user_input_phi = IntVar()
    user_input_phi.trace("w", lambda name, index, mode, var=user_input_phi: draw_angle_point_rep_canvas(rot_angles, [user_input_theta.get(), user_input_phi.get(), user_input_psi.get()]))
    slide_horizontal_1 = Scale(window, from_=0, to=360, orient=HORIZONTAL, variable=user_input_phi)  # this is not phi
    slide_horizontal_1.grid(column=3, row=10)

    slide_horizontal_3_text = Label(window, text="Angle around x (roll)", width="10")
    slide_horizontal_3_text.grid(column=2, row=11)
    user_input_psi = IntVar()
    user_input_psi.trace("w", lambda name, index, mode, var=user_input_psi: draw_angle_point_rep_canvas(rot_angles, [user_input_theta.get(), user_input_phi.get(), user_input_psi.get()]))
    slide_horizontal_3 = Scale(window, from_=0, to=360, orient=HORIZONTAL, variable=user_input_psi)  # this is not psi
    slide_horizontal_3.grid(column=3, row=11)

    octant_sel_bttn = Menubutton(window, text="Octant Selecter", width="40")
    octant_1 = BooleanVar()
    octant_2 = BooleanVar()
    octant_3 = BooleanVar()
    octant_4 = BooleanVar()
    octant_5 = BooleanVar()
    octant_6 = BooleanVar()
    octant_7 = BooleanVar()
    octant_8 = BooleanVar()

    octant_sel = np.array([True, True, True, True, True, True, True, True])
    Menu1 = Menu(octant_sel_bttn)
    Menu1.add_checkbutton(label="Octant 1", variable=octant_1)
    Menu1.add_checkbutton(label="Octant 2", variable=octant_2)
    Menu1.add_checkbutton(label="Octant 3", variable=octant_3)
    Menu1.add_checkbutton(label="Octant 4", variable=octant_4)
    Menu1.add_checkbutton(label="Octant 5", variable=octant_5)
    Menu1.add_checkbutton(label="Octant 6", variable=octant_6)
    Menu1.add_checkbutton(label="Octant 7", variable=octant_7)
    Menu1.add_checkbutton(label="Octant 8", variable=octant_8)
    octant_sel_bttn["menu"] = Menu1
    octant_sel_bttn.grid(sticky="W", columnspan="1", column=1, row=5)

    update_octant_bttn = Button(window, text='update octants',
                                command=lambda: update_rot_angles(n_sphere_points.get(),
                                                                  [user_input_theta.get(),
                                                                   user_input_phi.get(),
                                                                   user_input_psi.get()]))

    update_octant_bttn.grid(sticky="W", column=0, row=5)

    update_rot_angles(100, [0, 0, 0])
    # Let the window wait for any events
    window.mainloop()
