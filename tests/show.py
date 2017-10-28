from matplotlib import pyplot as plt


def showData():
	with open("data.txt") as f:
		data = [line.split() for line in f]
	x, y = zip(*data)
	plt.plot(x, y)
	plt.show()


if __name__ == "__main__":
	showData()
