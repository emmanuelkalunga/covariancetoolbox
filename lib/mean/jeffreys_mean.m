function JM = jeffreys_mean(COV)

E = mean_covariances(COV,'arithmetic');
H = mean_covariances(COV,'harmonic');

JM = riemann_geodesic(E,H,0.5);
