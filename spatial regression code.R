library(sp) #Calling library sp, You need to first install it
library(gstat) #Calling library gstat, You need to first install it
library(readxl)
library(expm) #Calling library expm, You need to first install it
year_1 = read_xlsx("~/spat1718.xlsx")
year_2 = read_xlsx("~/spat1819.xlsx")
year_3 = read_xlsx("~/spat1920.xlsx")

coordinates(year_2) = ~x+y
z2=year_2$deaths #My response variable is log(zinc) from dataset meuse1
#class(year_1)
v2 = variogram(z2~v1,year_2) #Obtaining variogram values with model of z vs lead
m2=fit.variogram(v2, vgm(c("Exp","Sph","Gau"))) #Fitting three types of variogram
m2 #Obtaining the best fitted variogram type
plot(v2, m2, main="Variogram fitting for the Diarrhoea Death Counts (Year 2018-2019)") #Plotting variogram
year_2$z1=year_2$x
year_2$z2=year_2$y


coordinates(year_1) = ~x+y
z1=year_1$deaths #My response variable is log(zinc) from dataset meuse1
#class(year_1)
v1 = variogram(z1~v1,year_1) #Obtaining variogram values with model of z vs lead
m1=fit.variogram(v1, vgm(c("Exp","Sph","Gau"))) #Fitting three types of variogram
m1 #Obtaining the best fitted variogram type
plot(v1, m1, main="Variogram fitting for the Diarrhoea Death Counts (Year 2017-18)") #Plotting variogram
year_1$z1=year_1$x
year_1$z2=year_1$y

coordinates(year_3) = ~x+y
z3=year_3$deaths #My response variable is log(zinc) from dataset meuse1
#class(year_1)
v3 = variogram(z3~v1,year_3) #Obtaining variogram values with model of z vs lead
m3=fit.variogram(v3, vgm(c("Exp","Sph","Gau"))) #Fitting three types of variogram
m3 #Obtaining the best fitted variogram type
plot(v3, m3, main="Variogram fitting for the Diarrhoea Death Counts (Year 2019-2020)") #Plotting variogram
year_3$z1=year_3$x
year_3$z2=year_3$y



meuse1$z1=meuse1$x #Defining variable to store latitude
meuse1$z2=meuse1$y #Defining variable to store longitude
coordinates(meuse1) = ~x+y #Combining both to obtain tuple of coordinates

z=log(meuse1$zinc) #My response variable is log(zinc) from dataset meuse1
v = variogram(z~1, meuse1) #Obtaining variogram values with model of z vs lead
plot(v, type="l") #Plotting variogram
m=fit.variogram(v, vgm(c("Exp","Sph","Gau","Mat"))) #Fitting three types of variogram
m #Obtaining the best fitted variogram type

# From the output we get that exponential type is best fitted variogram

c_0=0 #Nugget, We get this from output. 
c_e=2669.046 #Partial sill, We get this from output.
a_e=1.702381/3 #Range, We get this from output.

c_0=0 #Nugget, We get this from output. 
c_g=27715.39 #Partial sill, We get this from output.
a_g=1.918232/sqrt(3) #Range, We get this from output.

distance=matrix(nrow=36,ncol=36) #Creating a blank matrix to store distances
#For loop to obtain distances between every pair of centroids
for (i in 1:36) 
{
  for (j in 1:36)
  {
    distance[i,j]=sqrt((year_3$z1[i]-year_3$z1[j])^2+(year_3$z2[i]-year_3$z2[j])^2)
  }
}
for (i in 1:36)
{
  distance[i,i]=0
}

c_exp=matrix(nrow=36,ncol=36) #Creating a blank matrix to store variance-covariance matrix values
#For loop to obtain variance-covariance matrix values
for (i in 1:36)
{
  for (j in 1:36)
  {
    if(distance[i,j]==0)
    {
      c_exp[i,j]=c_0+c_e
    }
    else
    {
      c_exp[i,j]=c_0 + c_e * exp(-distance[i,j]/a_e)
    }
  }
}

c=matrix(nrow=36,ncol=36) #Creating a blank matrix to store variance-covariance matrix values
#For loop to obtain variance-covariance matrix values
for (i in 1:36)
{
  for (j in 1:36)
  {
    if(distance[i,j]==0)
    {
      c[i,j]=c_0+c_g
    }
    else
    {
      c[i,j]=c_g * exp(-(distance[i,j]/a_g)^2)
    }
  }
}





pow=solve(sqrtm(c)) #Obtaining C^{-1/2}

response_new=pow%*%z #Transforming response variable
explanatory_new=pow%*%year_1$v1 #Transforming Explanatory variable
model=lm(response_new~explanatory_new) #Fitting linear model on transformed variables
summary(model) #Obtaining summary of fitted linear model
plot(explanatory_new, response_new, main="2017-18 Rotavirus Vaccine Dose 1's response on Diarrhoea Death Counts", abline(lm(response_new~explanatory_new)), xlab="Vaccine Dose 1 Recipients (Transformed)", ylab="Diarrhoea Disease Counts (Transformed)")

#From the output you can obtain GLS estimates of beta's, sigma^2 and test the significance of beta's. You can also obtain R^2