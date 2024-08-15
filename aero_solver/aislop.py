import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import io

# Example function to create frames
def create_frames():
    frames = []
    t = np.linspace(0, 2 * np.pi, 100)
    for i in range(100):
        fig, ax = plt.subplots()
        y = np.sin(t + i * 0.1)
        ax.plot(t, y)
        ax.set_xlim(0, 2 * np.pi)
        ax.set_ylim(-1.5, 1.5)
        # Save the plot as an image in a buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        frames.append(buf)
        plt.close(fig)
    return frames

# Load frames as images
def load_frames(frame_buffers):
    images = []
    for buf in frame_buffers:
        image = plt.imread(buf)
        images.append(image)
    return images

# Example frames creation
frame_buffers = create_frames()
frames = load_frames(frame_buffers)

# Create animation
fig, ax = plt.subplots()



def update(frame):
    ax.clear()  # Clear the previous frame
    ax.imshow(frames[frame], animated=True)
    return []

ani = animation.FuncAnimation(fig, update, frames=len(frames), blit=True)

# Save the animation
ani.save('animation.mp4', writer='ffmpeg')

plt.show()