# Disroc variables
 Variables


- Calculation type: 

Ex.		!  Var   ! Value   
4           ! NB(12) | 1:Mechanics, 2:Hydraulic, 3:Thermal, 4:HM, 5:TM, 6:THM   
2           ! NB(15) | 1:Time Independent, 2:Transient   
1           ! NB(13) | 1: Plane, 2: Axismmetry   
0           ! var(1) | Symmetry Axis Position   
1           ! NB(14) | 1:Plane Strain 2: Stress   
0           ! var(2) | Generalized Strain or Stress   
1           ! NB(54) | 1:Permeable Matrix or 2:Dicrete Fracture Network   
0           ! NB(56) | 1: Gravity forces for hydraulic diffusion   
0.0098      ! var(16)| Vertical fluid pressure gradient due to gravity   
1           ! NB(72) | 1: User-defined process   


- Element type:

ntyp(n)     ! n element number.     
2           ! Bar element    
3           ! Triangular element   
4           ! Quatrilateral element   
5           ! Joint element   

# File results

## File numebr
|         |matrix mecha |matrix hydro | Fracture mecha | Fracture hydro | 
|--       |--           |--           |--              |--              |
|Element  |Fich101      | Fich103     | Fich201        |  Fich203       |
| Nodal   |Fich102      | Fich104     |                |                |
|Average  |             |             | Fich202        |                |



|        | matrix elem | matrix nodal | fracture elem  | fracture nodal |
| ----   |------------ | -----------  |----------------|----------------|
|Geometry|Fich301      |  Fiche302    |   Fiche 303    |                |



## File name
|         |matrix mecha     |matrix hydro      | Fracture mecha    | Fracture hydro | Geometry      |
|--       |--               |--                |--                 |--              | ---           |
|Element  |matrixMechaElem  | matrixHydroElem  | jointMecha        |  jointHydro    | matrixElemGeo |
| Nodal   |matrixMechaNodal | matrixHydroNodal |                   |                | nodeCoord     |
|Average  |                 |                  | jointMechaAverage |                |               |


|        | matrix elem | matrix nodal | fracture elem  | fracture nodal |
| ----   |------------ | -----------  |----------------|----------------|
|Geometry|matrixElemGeo|  nodeCoord   |   jointElemGeo |                |


## File structure

### Matrix results
- Fich101 matrixMechaElem     

| Temps du calcul |NoElem |EXX |EYY|EZZ|2EXY|Evol|SXX|SYY|SZZ|SXZ|meanStress|
|--               |    ---|--  |-- |-- |--  |--  |-- | --|-- |-- |---       |

- Fich102 matrixMechaNodal

|Temps du calcul|NoNode|Ux|Uy|
|---            |--    |--|--|

- Fich103 matrixHydroElem

|Temps du calcul|NoElem|Gausse pressure|vf|
|---            |--    |--|--|

- Fich104 matrixHydroNodal

|Temps du calcul|NoNode|Nodal pressure|
|---            |--    |--            |

### Joint results
- Fich201 jointMecha
  
|Temps du calcul|NoElem|Ut|Un|Tau|Sn|Utp|Unp|Damage variable|
|---            |   -- |--|--|-- |--|-- |-- |--             |


- Fich202 jointMecha
  
|Temps du calcul|NoElem|AUt|AUn|ATau|ASn|Average Damage|
|--             |---   |-- |-- |--  |-- |---           |

- Fich203 jointHydro

|Temps du calcul|NoElem|Gauss Pressure|Initial Pressure|
|--             |--    |--            |--              |

## Geometrical information

- Fich301
  
|MatrixElemNo|Surface|Node1|Node2|Node3|
|--          |--     |--   |--   |--   |

- Fiche302

|NodeNo|XNode|YNode|
|--    |--   |--   |

- Fich303

|JointElemNo|Length|Node1|Node2|Node3|Node4|
|--         |--    |--   |--   |--   |--   |




# Fracture Constitutive law, model 21510

- $beta$ controlling the damage varible $D$
  - A bigger $beta$ correspond to a more *ductile* fracture 
  - A smaller $beta$ correpond to a more *fragile* fracture
- $beta'$ controlling the friction angle $\phi$
  - A bigger $beta$ correspond to a more *fragile* fracture 
  - A smaller $beta$ correpond to a more *ductile* fracture

Example, testing cases, fracgile vs ductile fracture
All the other variables:

$k_t =2e8MPa, K_n=5e8MPa, e=10m,\sigma_R = 1.2MPa, C=5.8MPa, \phi=31,h_r=0.33$

|Case|$beta$|$beta'$|plasticity| $\sigma_n$|
|--- |---   |--     |--        |--         |
|1   | 1    |  1    | 0        |  -1MPa    |
|2   | 0.1  |  1    | 0        |  -1MPa    |
|3   | 1.8  |  1    | 0        |  -1MPa    |
|4   | 1    |  1    | 1        |  -1MPa    |
|5   | 0.1  |  1    | 1        |  -1MPa    |
|6   | 1.8  |  1    | 1        |  -1MPa    |
|7   | 1    |  1    | 0        |  -30MPa   |
|8   | 0.1  |  1    | 0        |  -30MPa   |
|9   | 1.8  |  1    | 0        |  -30MPa   |
|10  | 1    |  1    | 1        |  -30MPa   |
|11  | 0.1  |  1    | 1        |  -30MPa   |
|12  | 1.8  |  1    | 1        |  -30MPa   |

![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre1.png =400x) 
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre2.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre3.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre4.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre5.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre6.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre7.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre8.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre9.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre10.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre11.png =400x)
![pic](./Projects/Faultbehaviour/SingleFracture0/sensitivity/SingleFracutre12.png =400x)

# Fracture Constitutive law, Linear plastic with Mohr-Coulomb plasticity


