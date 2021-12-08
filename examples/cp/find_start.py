warnings.filterwarnings("ignore")

if __name__ == "__main__":
  data = np.loadtxt("LanlTi-sample.txt")
  print(data)

  print(np.shape(data))
  ifn = inter.Rbf(*(data[:,i] for i in range(data.shape[1])),
      function = 'inverse')
  
  space = [(0.0, 1.0),(0.0, 1.0),(0.0, 1.0),
                 (0.0, 1.0), (0.0, 1.0), (0.0, 1.0),
                 (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), 
                 (0.0, 1.0), (0.0, 1.0)]
  
  res = opt.minimize(lambda x: ifn(*x), [(i+j)/2 for i,j in space], 
      method = 'L-BFGS-B', bounds = space)
  if not res.success:
    warnings.warn("Surrogate model minimization did not succeed!")
  print(res.x)
  print(res.fun)