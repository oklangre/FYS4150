from matplotlib.pyplot import *

def make_plot(filename):
   infile = open(filename, 'r')
   infile.readline()       # skip first line
   x = []
   exact = []
   approx = []
   for line in infile:
      words = line.split()
      x_value = float(words[0])
      x.append(x_value)
      exact_value = float(words[2])
      exact.append(exact_value)
      approx_value = float(words[1])
      approx.append(approx_value)
   infile.close()
   return x, exact, approx

x10, u10, y10 = make_plot('general_thomas_10')
x100, u100, y100 = make_plot('general_thomas_100')
x1000, u1000, y1000 = make_plot('general_thomas_1000')
plot(x1000,u1000, x10, y10, x100, y100, x1000, y1000 )
legend(['Exact', 'Approximate n=10', 'Approximate n=100', 'Approximate n=1000'])
xlabel('x')
ylabel('u(x)') 
show()
print x
