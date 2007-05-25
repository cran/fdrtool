### gcmlcm.R  (2006-12-29)
###
###
###     Greatest Convex Minorant (GCM) and Least Concave Majorant (LCM) 
###
###
### Copyright 2006 Korbinian Strimmer 
###
###
### This file is part of the `fdrtool' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# find greatest convex minorant (gcm) or  
# least concave majorant (lcm)

gcmlcm = function(x, y, type=c("gcm", "lcm"))
{
  type=match.arg(type)
  
  if (is.unsorted(x))
    stop("The x values must be arranged in sorted order!")
  
  if (any(duplicated(x)))
    stop("No duplicated x values allowed!")

  ########
  
  dx = diff(x)
  dy = diff(y)

  rawslope = dy/dx

  if (type == "gcm") slope <- pvt.isoMean(rawslope, dx)
  if (type == "lcm") slope <- -pvt.isoMean(-rawslope, dx)

  # remove duplicate slopes
  keep = !duplicated(slope)
  x.knots = x[c(keep, TRUE)] # also keep last point
  dx.knots = diff(x.knots)
  slope.knots = slope[keep]

  y.knots = y[1]+c(0, cumsum(dx.knots*slope.knots))
  
  list(x.knots=x.knots,
       y.knots=y.knots,
       slope.knots=slope.knots)
}
