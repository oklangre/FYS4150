# program to find max value of relative error

def max_value(filename):
   infile = open(filename, 'r')
   infile.readline()
   error = []
   for line in infile:
      words = line.split()
      error_value = float(words[3])
      error.append(error_value)
   infile.close()
   return error

error_relative_list = max_value('thomasalgo_special_10000000')
error_relative = max(error_relative_list)
print error_relative, 10**error_relative
