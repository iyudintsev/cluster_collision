from matplotlib import pyplot as plt


def getInt(s):
	return int(s.split('=')[1])

def main():
	data = []
	with open("results.dat") as f:
		for line in f:
			parse_line = line.split(' ')
			data.append((getInt(parse_line[2]), 
						 getInt(parse_line[3])))
	data = list(zip(*data))
	plt.plot(data[0], data[1], marker='o')
	plt.show()


if __name__ == "__main__":
	main()
