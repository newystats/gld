# These are the expressions used in calculating the standard errors of the L-Moment estimates for the gld gpd

nu1 <- function(lambda){
  (2*lambda+1)*(2*lambda+3)*(2*lambda+5)*(2*lambda+7)*(lambda^2+2*lambda-5)^2
}

nu2 <- function(lambda){
  4*lambda^10 + 32*lambda^9 + 51*lambda^8 - 104*lambda^7 - 43*lambda^6 +
    735 * lambda^5 + 995 * lambda^4 + 6418 * lambda^3 + 22611*lambda^2 + 
    27911*lambda + 13966
}

nu3 <- function(lambda){
  (lambda+3)*(lambda^13 + 17*lambda^12 + 99*lambda^11 + 255*lambda^10 +
    667*lambda^9 + 3595*lambda^8 + 8745*lambda^7 - 17879*lambda^6 - 
    149808*lambda^5 - 312756*lambda^4 - 103684*lambda^3 + 
    432056*lambda^2 + 215458 * lambda - 367080)
}

nu4 <- function(lambda){
  (lambda+4)^2*(lambda^2-2*lambda -11)^2
}

nu5 <- function(lambda){
  lambda^8 + 8*lambda^7 + 12*lambda^6 - 84*lambda^5 - 405*lambda^4 -
    492*lambda^3 + 548*lambda^2 + 964*lambda - 840
}

nu6 <- function(lambda){
  2/(lambda*(lambda+1)*(lambda+2))
}

nu7 <- function(lambda){
    4*lambda^15 + 52*lambda^14 + 243*lambda^13 + 167*lambda^12 -
    2343*lambda^11 - 3244*lambda^10 + 27697*lambda^9 +
    74517*lambda^8 + 8739*lambda^7 + 244921*lambda^6 +
    1532136*lambda^5 + 1771024*lambda^4 - 923048*lambda^3 -
    1541641*lambda^2 + 686268*lambda + 550620
}

nu8 <- function(lambda){
  4*lambda^9 + 32*lambda^8 +35*lambda^7 - 216*lambda^6 - 159*lambda^5 +
    1053*lambda^4 + 350*lambda^3 - 1451*lambda^2 + 148*lambda + 420
}

nu9 <- function(lambda){
  8*lambda^10 + 124*lambda^9 + 822*lambda^8 + 3133*lambda^7 + 
    7894*lambda^6 + 14519 *lambda^5 + 22236*lambda^4 + 35026*lambda^3 + 
    53894*lambda^2 + 54316*lambda + 22696
}

nu10 <- function(lambda){
  4*lambda^14 + 84*lambda^13 + 757*lambda^12 + 3614*lambda^11 + 
    7369*lambda^10 - 18988*lambda^9 - 199465*lambda^8 - 
    747082*lambda^7 - 1606417*lambda^6 - 1893644*lambda^5 - 
    580072*lambda^4 + 1264048*lambda^3 + 1118160*lambda^2 -
    396656*lambda - 567840
}

nu11 <- function(lambda){
    8*lambda^15 + 172*lambda^14 +1742*lambda^13 + 10633*lambda^12 +
    39856*lambda^11 + 79701*lambda^10 + 31792*lambda^9 -
    125933*lambda^8 + 172824*lambda^7 + 1449179*lambda^6 +
    2083778*lambda^5 + 200988*lambda^4 - 1588440*lambda^3 -
    552332*lambda^2 + 624208*lambda +283920
}

nu12 <- function(lambda){
  18*(lambda+1)*(lambda+3)*(2*lambda^2+3*lambda+2)
}

nu13 <- function(lambda){
  2*lambda^10 + 29*lambda^9 + 149*lambda^8 + 253*lambda^7 - 575*lambda^6 - 
  3857*lambda^5 - 8701*lambda^4 - 6229*lambda^3 + 11573*lambda^2 + 14412*lambda - 15120
}  

nu14 <- function(lambda){
    29*lambda^7 + 206*lambda^6 - 119*lambda^5 - 3822*lambda^4 - 7189*lambda^3 +
      5506*lambda^2 + 13279*lambda -10290
}
  
nu15 <- function(lambda){
    4*lambda^11 + 28*lambda^10 - 49*lambda^9 - 699*lambda^8 - 1280*lambda^7 -
      1601*lambda^6 - 11409*lambda^5 - 22185*lambda^4 + 5456*lambda^3 +
      23485*lambda^2 - 6654*lambda - 7560
}
  
nu16 <- function(lambda){
    8*lambda^13 + 44*lambda^12 - 282*lambda^11 - 1699*lambda^10 + 4418*lambda^9 +
      32076*lambda^8 - 1094*lambda^7 - 164198*lambda^6 + 24618*lambda^5 +
      557596*lambda^4 - 7256*lambda^3 - 643599*lambda^2 + 100548*lambda +
      185220
}

nu17 <- function(lambda){
    4*lambda^5 + 4*lambda^4 - 17*lambda^3 + 43*lambda^2 - 14*lambda +16
}
