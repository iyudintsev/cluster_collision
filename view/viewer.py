import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d


class DataViewer(object):
	def __init__(self, L, pause):
		self.L = L
		self.file = open("data.txt")
		self.pause = pause

	def __del__(self):
		self.file.close()

	def get(self):
		data = []
		while True:
			line = self.file.readline()
			if line.strip() == "next":
				break
			if not line:
				break
			data.append([float(num) for num in line.split()])

		return data

	def run(self):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		
		data = self.get()
		while data:
			plt.cla()
			x, y, z = list(zip(*data))
			ax.scatter(x, y, z, c='m', marker='o')
			
			ax.set_xlim(0, self.L)
			ax.set_ylim(0, self.L)
			ax.set_zlim(0, self.L)

			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_zlabel('z')
			ax.set_title("MD")

			plt.pause(self.pause)
			data = self.get()


def main():
	data_viewer = DataViewer(20, 0.01)
	data_viewer.run()


if __name__ == "__main__":
	main()
