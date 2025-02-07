{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "This is the Python code of Example 3.13 in the book:\n",
    "\n",
    "FIBER-REINFORCED COMPOSITES\n",
    "Materials, Manufacturing, and Design\n",
    "\n",
    "by: P.K. Mallick (2008) by Taylor & Francis Group, LLC\n",
    "\n",
    "This code:\n",
    "(C) Joris Remmers (2013-2019)\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Example 3.13"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Calculate lamina stresses at the midplane of each lamina in the [+45/-45] laminate in Example 3.7 due to $N_{xx} = 100 kN/m$."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Solution\n",
    "\n",
    "Import the correct functions from the composite module and required mathermatical operators:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "from composite    import TransverseIsotropic,mixMaterials,Laminate,stressTransformation\n",
    "from numpy        import array,dot,zeros\n",
    "from numpy.linalg import inv"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "As demonstrated in Example 3.6 and 3.7, the ply properties and laminate can be obtained as: "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "  Laminate properties\n",
      "  -----------------------------------------------------------\n",
      "  layer   thick orient.  material\n",
      "  -----------------------------------------------------------\n",
      "      0   0.006    -45   composite\n",
      "      1   0.006     45   composite\n",
      "\n"
     ]
    }
   ],
   "source": [
    "carbon = TransverseIsotropic( 220e9,0.2,91.7e9)\n",
    "epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)\n",
    "\n",
    "compmat = mixMaterials( carbon , epoxy , 0.6 )\n",
    "\n",
    "lam = Laminate()\n",
    "\n",
    "lam.addMaterial( 'composite' , compmat )\n",
    "\n",
    "lam.addLayer( 'composite' , -45 , 6e-3 )\n",
    "lam.addLayer( 'composite' , 45 , 6e-3 )\n",
    "\n",
    "print(lam)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "To be able to calculate the midplain strains and curvature, the inverse laminate stiffness matrices $[A_1]$, $[B_1]$, $[C_1]$, and $[D_1]$ are required. First, the $[A]$, $[B]$ and $[D]$ matrixes are computed, and then the corresponding inverse matrice."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "The matrix A1:\n",
      " [[ 7.73727620e-09 -5.06668026e-09  0.00000000e+00]\n",
      " [-5.06668026e-09  7.73727620e-09  0.00000000e+00]\n",
      " [ 0.00000000e+00  0.00000000e+00  5.69566724e-09]] \n",
      "\n",
      "The matrix B1:\n",
      " [[-0.0000000e+00 -0.0000000e+00 -6.0459327e-07]\n",
      " [-0.0000000e+00 -0.0000000e+00 -6.0459327e-07]\n",
      " [-6.0459327e-07 -6.0459327e-07 -0.0000000e+00]] \n",
      "\n",
      "The matrix C1:\n",
      " [[-0.0000000e+00 -0.0000000e+00 -6.0459327e-07]\n",
      " [-0.0000000e+00 -0.0000000e+00 -6.0459327e-07]\n",
      " [-6.0459327e-07 -6.0459327e-07 -0.0000000e+00]] \n",
      "\n",
      "The matrix D1:\n",
      " [[ 0.00064477 -0.00042222  0.        ]\n",
      " [-0.00042222  0.00064477  0.        ]\n",
      " [ 0.          0.          0.00047464]] \n",
      "\n"
     ]
    }
   ],
   "source": [
    "A = lam.getA()\n",
    "B = lam.getB()\n",
    "D = lam.getD()\n",
    "\n",
    "A1,B1,C1,D1 = lam.getInverseMatrices()\n",
    "\n",
    "print(\"The matrix A1:\\n\",A1,\"\\n\")\n",
    "print(\"The matrix B1:\\n\",B1,\"\\n\")\n",
    "print(\"The matrix C1:\\n\",C1,\"\\n\")\n",
    "print(\"The matrix D1:\\n\",D1,\"\\n\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Now the midplain strains and curvature can be calculated: "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "The midplane strains are: [ 0.00077373 -0.00050667  0.        ]\n",
      "The curvatures are: [ 0.          0.         -0.06045933]\n"
     ]
    }
   ],
   "source": [
    "N = zeros(3)\n",
    "N[0] = 1e5\n",
    "\n",
    "eps0 = dot(A1,N)\n",
    "\n",
    "print(\"The midplane strains are:\",eps0)\n",
    "\n",
    "kappa = dot(C1,N)\n",
    "\n",
    "print(\"The curvatures are:\",kappa)\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Using these midplain strains and curvature, the stresses and strains in the individual layers can be calculated as:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "\n",
      "The strains in the bottom layer : \n",
      " [ 0.00077373 -0.00050667  0.00018138]\n",
      "\n",
      "The strains in the top layer    : \n",
      " [ 0.00077373 -0.00050667 -0.00018138]\n",
      "\n",
      "The stresses in the bottom layer : \n",
      " [ 8.33333333e+06  9.31322575e-10 -2.08995545e+06]\n",
      "\n",
      "The stresses in the top layer    : \n",
      " [8.33333333e+06 9.31322575e-10 2.08995545e+06]\n"
     ]
    }
   ],
   "source": [
    "epsb = eps0-3e-3*kappa\n",
    "epst = eps0+3e-3*kappa\n",
    "\n",
    "print (\"\\nThe strains in the bottom layer : \\n\",epsb)\n",
    "print (\"\\nThe strains in the top layer    : \\n\",epst)\n",
    "\n",
    "stressb = dot(lam.getQbar(0),epsb)\n",
    "stresst = dot(lam.getQbar(1),epst)\n",
    "\n",
    "print (\"\\nThe stresses in the bottom layer : \\n\",stressb)\n",
    "print (\"\\nThe stresses in the top layer    : \\n\",stresst)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The stresses can be expressed in the material frame of reference (12) by transformation:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "-0.7853981633974483\n",
      "\n",
      "The stress in the bottom layer : \n",
      " [6256622.12152551 2076711.21180782 4166666.66666666]\n",
      "0.7853981633974483\n",
      "\n",
      "The stress in the top layer    : \n",
      " [ 6256622.12152551  2076711.21180782 -4166666.66666666]\n"
     ]
    }
   ],
   "source": [
    "print(\"\\nThe stress in the bottom layer : \\n\",stressTransformation(stressb,-45.))\n",
    "print(\"\\nThe stress in the top layer    : \\n\",stressTransformation(stresst,45.))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
