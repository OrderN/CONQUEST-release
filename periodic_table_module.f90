module periodic_table

  use datatypes

  implicit none
  
  integer, parameter  :: n_species=94

  real(double), parameter ::  atomic_mass(n_species) = (/ &
       1.01_double,    4.00_double,  6.94_double,  9.01_double, 10.81_double, 12.01_double, &
       14.01_double,  16.00_double, 19.00_double, 20.18_double, 22.99_double, 24.31_double, &
       26.98_double,  28.09_double, 30.97_double, 32.07_double, 35.45_double, 39.95_double, &
       39.10_double,  40.08_double, 44.96_double, 47.88_double, 50.94_double, 52.00_double, &
       54.94_double,  55.85_double, 58.93_double, 58.69_double, 63.55_double, 65.39_double, &
       69.72_double,  72.61_double, 74.92_double, 78.96_double, 79.90_double, 83.80_double, & 
       85.47_double,  87.62_double, 88.91_double, 91.22_double, 92.91_double, 95.94_double, & 
       98.91_double, 101.07_double,102.91_double,106.42_double,107.87_double,112.41_double, &
       114.82_double,118.71_double,121.75_double,127.60_double,126.90_double,131.29_double, &
       132.91_double,137.33_double,138.91_double,140.12_double,140.91_double,144.24_double, &
       146.92_double,150.36_double,151.97_double,157.25_double,158.93_double,162.50_double, &
       164.93_double,167.26_double,168.93_double,173.04_double,174.97_double,178.49_double, &
       180.95_double,183.85_double,186.21_double,190.20_double,192.22_double,195.08_double, &
       196.97_double,200.59_double,204.38_double,207.20_double,208.98_double,208.98_double, &
       209.99_double,222.02_double,223.02_double,226.03_double,227.03_double,232.04_double, &
       231.04_double,238.03_double,237.05_double,244.06_double/)
  
end module periodic_table
