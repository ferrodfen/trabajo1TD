###### Función para aplicar el Método de Hurwicz que nos permita obtener los valores de alpha para los que cambian las alternativas óptimas

## Funciones que vamos a utilizar del fichero "teoriadecision_funciones_incertidumbre.R" (he copiado y pegado)

criterio.Hurwicz = function(tablaX,alfa=0.3,favorable=TRUE) {
  # alfa es un escalar entre 0 y 1 lo obtiene para ese único valor
  X = tablaX;
  if (favorable) {
    Altmin = apply(X,MARGIN=1,min);
    Altmax= apply(X,MARGIN=1,max);
    AltH = alfa * Altmax + (1-alfa) * Altmin 
    Hurwicz = max(AltH)
    Alt_Hurwicz = which.max.general(AltH)
    metodo = 'favorable';
  } else {
    Altmin = apply(X,MARGIN=1,min);
    Altmax= apply(X,MARGIN=1,max);
    AltH = (1-alfa) * Altmax + alfa * Altmin 
    Hurwicz = min(AltH)
    Alt_Hurwicz = which.min.general(AltH)
    metodo = 'desfavorable';
  }
  resultados = list();
  resultados$criterio = 'Hurwicz';
  resultados$alfa = alfa;
  resultados$metodo = metodo;
  resultados$tablaX = tablaX;
  resultados$ValorAlternativas = AltH;
  resultados$ValorOptimo = Hurwicz;
  resultados$AlternativaOptima = Alt_Hurwicz;
  
  return(resultados);
  
}



dibuja.criterio.Hurwicz = function(tablaX,favorable=TRUE){
  X = tablaX;
  Altmin = apply(X,MARGIN=1,min);
  Altmax = apply(X,MARGIN=1,max);
  valfa = seq(from=0,to=1,by=0.05);
  vHurwicz = rep(0,length(valfa));
  Alt_vHurwicz = rep(0,length(valfa));
  for (i in 1:length(valfa)) {
    alfab = valfa[i];  
    if (favorable) { 
      vAltH = alfab * Altmax + (1-alfab) * Altmin; 
      vHurwicz[i] = max(vAltH)
    } else {
      vAltH = alfab * Altmin + (1-alfab) * Altmax; 
      vHurwicz[i] = min(vAltH)      
      }
    }
  x0=0;x1=1;
  y0 = min(Altmin);
  y1 = max(Altmax);
  rg = y1-y0;
  y0=y0-0.1*rg;y1=y1+0.1*rg;
  plot(c(x0,x1), c(y0,y1), type = "n", xlab = "alpha", ylab = "Criterio Hurwicz"); 
  nn = length(Altmin);
  colores = rainbow(nn);
  abline(v=0);
  abline(v=1);
  if (favorable) {
    for (i in 1:nn) {
      aa = Altmin[i];
      bb = (Altmax[i] - Altmin[i]);
      abline(a=aa,b=bb,col=colores[i]);
    }
  } else {
    for (i in 1:nn) {
      aa = Altmax[i];
      bb = (Altmin[i] - Altmax[i]);
      abline(a=aa,b=bb,col=colores[i]);
    }        
  }
  lines(valfa,vHurwicz,col=rainbow(nn+1)[nn+1],lty=3,lwd=3)
  if (favorable) {
    legend("bottomright",legend=rownames(X),fill=colores,inset=0.05)
    title("Criterio de Hurwicz (favorable - línea discontinua)")
  } else {
    legend("topright",legend=rownames(X),fill=colores,inset=0.05)
    title("Criterio de Hurwicz (desfavorable - línea discontinua)")    
  }

}



# FUNCION : esta funcion nos da los avlores de alfa para los que las alternativas cambian 

# Entrada: Tabla, favorable (T/F)
# Salida: Intervalo -> Alternativa (óptima para ese intervalo de alfa)


Hurwicz.intervalos = function(tabla, favorable = TRUE){
  X = tabla # renombramos la tabla
  Altmin = apply(X,MARGIN=1,min)      # vector de minimos (por filas)
  Altmax = apply(X,MARGIN=1,max)      # vector de maximos (por filas)
  valfa = seq(from=0,to=1,by=0.05)    # vector de valores para alfa (entre 0 y 1)
  Alt_opt = rep(0,length(valfa))      # creamos el vector de decisiones (por el criterio de Hurwicz) para cada valor de alfa
  alfaCorte=c()                       # vector que contiene los valores de alfa donde cambian las decisiones
  for (i in 1:length(valfa)) {
    Alt_opt[i] = criterio.Hurwicz(X, alfa = valfa[i], favorable)$AlternativaOptima # obtenemos las alternativas para cada alfa
    Alt=c() # Este va a ser el Vector de las alternativas optimas para todos los alfa
    for (i in 1:length(Alt_opt)) {
      valrepetidos = duplicated(Alt_opt) # Vector de TRUE/FALSE donde los FALSE son los elementos que se repiten 
      if (isFALSE(valrepetidos[i])){
        Alt = c(Alt,Alt_opt[i]) # Si es falso (si el valor se repite) lo almacenamos en el vector Alt
      }
    }
  }
  # Teniendo el vector de alternativas (Alt) buscamos los puntos de corte de las rectas asociadas a cada alternativa (beneficios)
  # Por ejemplo, la recta que sale de la alternativa a1 y a2 seria: 
  #
  #               a1Max *alfa +(1-alfa)*a1Min = a2Max *alfa +(1-alfa)*a2Min
  # 
  # Pasando todo a un  miembro e igualando a 0 nos queda: 
  #
  #               alfa * (a1Max- a2Max - a1Min + a2Min) + a1Min -a2Min = 0
  # 
  # Buscamos ahora los valores de alfa para los que se cortan las rectas asociadas a cada decision
  for (i in 1:(length(Alt)-1)){
    imax = as.numeric(Altmax[Alt[i]])      # maximo asociado a la decision i del vector Alt
    imax1 = as.numeric(Altmax[Alt[i+1]])   # maximo asociado a la decision i+1 del vector Alt
    imin = as.numeric(Altmin[Alt[i]])      # minimo asociado a la decision i del vector Alt
    imin1 = as.numeric(Altmin[Alt[i+1]])   # minimo asociado a la decision i+1 del vector Alt
    if (favorable){
      pCorte = function(alfa) {alfa * (imax-imax1-imin+imin1)+imin-imin1}
      alfaC = uniroot(pCorte, interval = c(0,1))$root[[1]] # Buscamos los 0 para cada funcion
      alfaCorte[i] = alfaC  # Almacenamos los valores de alfa para los que las rectas se cortan en alfaCorte
    } else {
      # Para el caso de costes (alternativas a1 y a2):
      # 
      #               a1Max *(1-alfa) +alfa*a1Min = a2Max *(1-alfa) +alfa*a2Min
      # 
      # Pasando todo a un  miembro e igualando a 0 nos queda: 
      #
      #               alfa * (a1Min- a2Min - a1Max + a2Max) + a1Max -a2Max = 0
      #
      pCorte = function(alfa) {alfa * (imin-imin1-imax+imax1)+imax-imax1}
      alfaC = uniroot(pCorte, interval = c(0,1))$root[[1]] # Buscamos los 0 para cada funcion
      alfaCorte[i] = alfaC  # Almacenamos los valores de alfa para los que las rectas se cortan en alfaCorte
    }
    
  }
  return(alfaCorte)
}


# Ejemplos para probar 
tb07 = crea.tablaX(c(45,45,540,-180,900,-900), numalternativas = 3, numestados = 2)
tb07

tb05 = crea.tablaX(c(125,120,156,60,130,80), numalternativas = 3, numestados = 2)
tb05
Hurwicz.intervalos(tb07, favorable)
Hurwicz.intervalos(tb05, favorable=F)
