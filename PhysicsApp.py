from tkinter import *
from tkinter import ttk
import tkinter as tk
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class PhysicsApp(Tk):
    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        self.title("Spring Motion")

        frame = Frame(self)
        frame.pack(side=TOP, padx=5, pady=5)

        # Create input fields and labels
        self.mass_label = tk.Label(frame, text="Mass (kg):")
        self.mass_label.grid(row=0, column=0, padx=5, pady=5)
        self.mass_entry = tk.Entry(frame)
        self.mass_entry.grid(row=0, column=1, padx=5, pady=5)

        self.spring_constant_label = tk.Label(frame, text="Spring constant (N/m):")
        self.spring_constant_label.grid(row=1, column=0, padx=5, pady=5)
        self.spring_constant_entry = tk.Entry(frame)
        self.spring_constant_entry.grid(row=1, column=1, padx=5, pady=5)

        self.initial_displacement_label = tk.Label(frame, text="Initial displacement (m):")
        self.initial_displacement_label.grid(row=2, column=0, padx=5, pady=5)
        self.initial_displacement_entry = tk.Entry(frame)
        self.initial_displacement_entry.grid(row=2, column=1, padx=5, pady=5)

        # Create button to calculate and display results
        self.calculate_button = tk.Button(frame, text="Calculate", command=self.calculate_spring_motion)
        self.calculate_button.grid(row=3, column=0, padx=5, pady=5)

        # Create labels to display results
        self.amplitude_label = tk.Label(frame, text="Amplitude:")
        self.amplitude_label.grid(row=4, column=0, padx=5, pady=5)
        self.amplitude_value = tk.Label(frame, text="")
        self.amplitude_value.grid(row=4, column=1, padx=5, pady=5)

        # self.period_label = tk.Label(frame, text="Period:")
        # self.period_label.grid(row=5, column=0, padx=5, pady=5)
        # self.period_value = tk.Label(frame, text="")
        # self.period_value.grid(row=5, column=1, padx=5, pady=5)

        self.equation_label = tk.Label(frame, text="Equation:")
        self.equation_label.grid(row=6, column=0, padx=5, pady=5)
        self.equation_value = tk.Label(frame, text="")
        self.equation_value.grid(row=6, column=1, padx=5, pady=5)

    def calculate_spring_motion(self):
        mass = float(self.mass_entry.get())
        spring_constant = float(self.spring_constant_entry.get())
        initial_displacement = float(self.initial_displacement_entry.get())

        # Calculate results
        omega = math.sqrt(spring_constant / mass)
        period = 2 * math.pi / omega
        amplitude = abs(initial_displacement)

        # Construct equation string
        if initial_displacement > 0:
            equation = f"f(x) = {amplitude:.2f} * sin({omega:.2f}t)"
        elif initial_displacement < 0:
            equation = f"f(x) = -{amplitude:.2f} * sin({omega:.2f}t)"
        else:
            equation = f"f(x) = 0"

        # Update labels with results
        self.amplitude_value.config(text=f"{amplitude:.2f}");
        # Update the period value label
        self.period_value.config(text=f"{period:.2f}");

        # Create a figure for the plot
        fig = Figure(figsize=(5, 4), dpi=100);
        ax = fig.add_subplot(111);

        # Generate the x values for the plot
        x_vals = np.linspace(0, 2 * np.pi, 200);

        # Generate the y values for the plot
        y_vals = amplitude * np.sin(omega * x_vals);

        # Add the plot to the figure
        ax.plot(x_vals, y_vals);

        # Add a title and labels to the plot
        ax.set_title("Harmonic Function");
        ax.set_xlabel("x");
        ax.set_ylabel("f(x)");

        # Create a canvas for the figure and add it to the GUI
        # canvas = FigureCanvasTkAgg(fig, frame=self);
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw();
        canvas.get_tk_widget().grid(sticky="nsew")
            #row=7, columnspan=2)

        # Period value label
        self.period_label = tk.Label(master=frame1, text="Period (s):")
        self.period_label.grid(row=2, column=0, sticky="w")
        self.period_value = tk.Label(master=frame1, text="0.00")
        self.period_value.grid(row=2, column=1, sticky="w")

        # Equation of motion label
        self.equation_label = tk.Label(master=frame, text="Equation of motion:")
        self.equation_label.grid(row=3, column=0, sticky="w")
        self.equation_value = tk.Label(master=frame1, text="f(x) = 0")
        self.equation_value.grid(row=3, column=1, sticky="w")

        # Graph canvas
        self.graph_canvas = tk.Canvas(master=frame2, width=400, height=400)
        self.graph_canvas.pack()

        # Plot button
        self.plot_button = tk.Button(master=frame3, text="Plot", command=self.plot_graph)
        self.plot_button.pack()

    def plot_graph(self):
        # Clear the canvas
        self.graph_canvas.delete("all")

        # Get the values from the text fields
        mass = float(self.mass_entry.get())
        k = float(self.spring_constant_entry.get())
        x0 = float(self.initial_position_entry.get())

        # Calculate the necessary values
        amplitude, period, equation = spring_motion(mass, k, x0)
        x_values = [i / 100 for i in range(-500, 501)]
        y_values = [harmonic_motion(x, amplitude, period, equation) for x in x_values]

        # Scale the x and y values to fit the canvas
        x_scale = 150
        y_scale = 150
        x_offset = 200
        y_offset = 200
        scaled_x_values = [x * x_scale + x_offset for x in x_values]
        scaled_y_values = [-y * y_scale + y_offset for y in y_values]

        # Draw the x and y axes
        self.graph_canvas.create_line(x_offset, 0, x_offset, 400)
        self.graph_canvas.create_line(0, y_offset, 400, y_offset)

        # Draw the function
        for i in range(len(scaled_x_values) - 1):
            x1 = scaled_x_values[i]
            y1 = scaled_y_values[i]
            x2 = scaled_x_values[i + 1]
            y2 = scaled_y_values[i + 1]
            self.graph_canvas.create_line(x1, y1, x2, y2, fill="blue")

        # Calculate the points to plot the harmonic function
        points = []
        for i in range(0, 1000):
            x = i * (xmax - xmin) / 1000 + xmin
            y = spring_motion.harmonic_motion(x, amplitude, period, phase_shift)
            points.append((x, y))

        # Plot the points on the canvas
        for i in range(len(points) - 1):
            x1, y1 = points[i]
            x2, y2 = points[i + 1]
            # Scale the points to fit the canvas
            x1 = (x1 - xmin) / (xmax - xmin) * plot_width + plot_padding
            y1 = plot_height - (y1 - ymin) / (ymax - ymin) * plot_height + plot_padding
            x2 = (x2 - xmin) / (xmax - xmin) * plot_width + plot_padding
            y2 = plot_height - (y2 - ymin) / (ymax - ymin) * plot_height + plot_padding
            self.graph_canvas.create_line(x1, y1, x2, y2, fill="blue")

        self.graph_canvas.create_text(self.canvas_width / 2, 20, text="Harmonic Function", font=("Arial", 16))

