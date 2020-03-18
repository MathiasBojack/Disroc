# Material infromation

    Material.
            project_path
            project_name
            total_number
            type{i}.
                    name
                    nature
                    mecha.
                        modelnum
                        numPara
                        Para(1:numPara)
                    hydro.
                        modelnum
                        numPara
                        Para(1:numPara)
                    couplingPar(1:3)

# Calculation parameter
    Parameter.
        problem.
                physics
                time
                type
                axesymmetry
                planetype
                generalized
                hyro.
                    matrix
                    gravity.active
                    gravity.value
                user

        specpara.
                staging
                stepnum
                boundaryforce
        
        load.
            resumption.active
            maxratio
            volumeforce.
                        active
                        gx
                        gy
    
            resumption.
                        stepnum    
                        stepactive 

        calpara.
                loadincrement
                itermax
                tolerance.
                        criteria
                        convergence
                        displacement
                time.
                    start
                    end
                    increment

# readDatFile

    prob_info.
              proj_name
              proj_path   containing '.gid'
              coordinate
              connectivity.
                           matrix
                           joint
              numNode
              numElem
              numJointElem
              numMat


## structure of Disroc file
|Start| End| Content|
|---  |--- |---     |
|1    |    | Disroc5|    
|2    |    |General information|    
|3    |    |Coordinates|    
|4    |numNode+3|Nodal coordinates|
|numNode+4||Connectivity|
|numNode+5|numNode+4+numElem|Connectivity of elements|
|numNode+5+numElem||Materials|
|numNode+6+numElem|numNode+5+numElem+numMat| Material info|
|numNode+6+numElem+numMat||Boundary|
|||Boundaryconditions|
|||End Data|





