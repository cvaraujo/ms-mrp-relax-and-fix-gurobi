import os
import subprocess

folders = ['50', '75', '100']
algorithms = ['model']
commands = []

delay = []
jitter = []
variation = []

for alg in algorithms:
	os.mkdir("results_" + alg)
os.mkdir("preprocessing")

for f in folders:
	for i in range(10, 110, 10):
		path = "Washington-" + f + "/" + str(i)
		instance = os.listdir(path)
		for inst in instance:
			if (inst[0] == 'w' or inst[0] == '1'):
				for alg in algorithms:
					commands.append("./MaxService " + path + "/" + inst + " " + path + "/param-" + inst + " results_" + alg + "/result_" + inst + " 1 1 " + alg + " 3600 preprocessing/" + inst)

f = '200'
for i in range(125, 375, 25):
	path = "Washington-" + f + "/" + str(i)
	instance = os.listdir(path)
	for inst in instance:
		if (inst[0] == 'w' or inst[0] == '1'):
			for alg in algorithms:
				commands.append("./MaxService " + path + "/" + inst + " " + path + "/param-" + inst + " results_" + alg + "/result_" + inst + " 1 1 " + alg + " 3600 preprocessing/" + inst)

for c in commands:
   print(c)
   p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
   msg, err = p.communicate()
   if msg:
   	print(msg)
   print("OK!!")
