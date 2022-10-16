import { readFileSync } from 'fs'
import { inspect } from 'util'
import 'util'
import v from 'vector-math'
const { createVectorObj, crossProduct, dotProduct, subVector, unitVector } = v



const PDB = readFileSync('./assets/1ul1.cif', 'utf-8').split('\n')
const ATOMS = PDB.filter((line) => line.startsWith('ATOM')).map(array => array.split(' ').filter(element => element !== ''))
const atomLabels = PDB.filter((line) => line.startsWith('_atom_site')).filter((line) => line.startsWith('_atom_sites.fract_transf') === false && line.startsWith('_atom_sites.entry_id') === false)
const customLabels = ['', 'id', 'atom', 'atom_type', '', 'residue', 'chain', 'entity_index', 'residue_index', '', 'x', 'y', 'z', 'occupancy', 'isotropic_temperature_factor', 'formal_charge', 'author_residue_index', 'author_residue', 'author_chain', 'author_atom_type', '' ]

let object = {}
const newAtoms = ATOMS.map((array) => {
    let atom_info = {}
    let author_entries = {}
    array.forEach((entry, index) => {
        const label = customLabels[index]
        if(label !== '') {
            if(label.includes('author')) {
                if( label.includes('index')) {
                    const int = parseInt(entry)
                    author_entries = {...author_entries, [label.replace('author_', '')]: int}
                }
                else {
                    author_entries = {...author_entries, [label.replace('author_', '')]: entry}
                }
                atom_info = {...atom_info, author_entries}
                return
            }
            if(label === 'id' || label.includes('index') || label === 'formal_charge') {
                const int = parseInt(entry)
                atom_info = ({...atom_info, [label]: int})
                return
            }
            if(label === 'isotropic_temperature_factor' || label === 'x' || label === 'y' || label === 'z' || label === 'occupancy') {
                const float = parseFloat(entry)
                atom_info = ({...atom_info, [label]: float})
                return
            }
            else {
                atom_info = ({...atom_info, [label]: entry})
            }
        }
    })
    return (
        atom_info
    )
})
const entityStart = PDB.indexOf('_entity.details ')
const chainData = []

for(let i=entityStart; i <= 1000; i++) {
    if(PDB[i] === '# ') {
        break
    } else {
        chainData.push(PDB[i])
    }
}

chainData.shift()

const chainInfo = chainData.map((line) => replaceSpaces(line, 0)).map((line) => replaceQuotes(line, 0)).map((line) => line.split(' ')).map((array) => array.filter((element) => element !== ''))

function replaceSpaces(line, count) {
    if(count === line.length-1) {
        return (
            line
        )
    }
    const newLine = line.replace(/(?<='\w+)\s/, '_')
    count++
    return replaceSpaces(newLine, count)
}

function replaceQuotes(line, count) {
    if(count >= 2) {
        return (
            line
        )
    }
    const newLine = line.replace(/['"`]/, '')
    count++
    return replaceQuotes(newLine, count)
}


const chainArray = chainInfo.map((entry) => {
    let chainObj = {}
    entry.forEach((entry, index) => {
        switch (index) {
            case 1:
                chainObj = {...chainObj, type: entry }
                break;
            case 3:
                chainObj = {...chainObj, name: entry }
                break;
            case 4:
                chainObj = {...chainObj, molecular_weight: parseFloat(entry) }
                break;
            case 5:
                chainObj = {...chainObj, quantity: parseInt(entry) }
                break;
        
            default:
                break;
        }
    })
    return chainObj
})


object = {...object, chain_info: chainArray}

let chainCount = 0
object.chain_info.filter((entry) => entry.type === 'polymer').forEach((entry) => chainCount += entry.quantity)

const chains = []
for(let i=0; i < chainCount; i++) {
    chains.push([])
}

const chainLabels = newAtoms.map((atom) => atom.chain).filter((chainID, index, array) => chainID !== array[index+1])
const chainAtoms = chains.map((chain, index) => newAtoms.filter((atom) => atom.chain === chainLabels[index]))


// seq with holes begin
const sequenceStart = PDB.indexOf('_entity_poly_seq.hetero ')
const seq = []

for(let i=sequenceStart; i <= 100000; i++) {
    if(PDB[i] === '# ') {
        break
    } else {
        seq.push(PDB[i])
    }
}

// console.log(chainAtoms[0])
const polymers = chainInfo.filter((chain) => chain[1] === 'polymer')
const sequences = polymers.map((chain, index) => seq.filter((line) => line.startsWith(index+1)).map((line) => line.split(' ').filter((entry) => entry !== '')[2]))

// console.log(sequences[0])
const seqChains = chainLabels.map((chain, index) => sequences[chainAtoms[index][0].entity_index-1])

const seqResidues = chainAtoms.map((chain, index) => {
        const residueArray = [[]]
        let n = 0
        chain.forEach((atom, index, array) => {
            if(index !== array.length-1) {
                if(atom.residue_index === array[index+1].residue_index) {
                    residueArray[n].push(atom)
                } else {
                    residueArray.push([])
                    n++
                }
            } else { residueArray[n].push(atom) }
        })
    return (
        residueArray
    )
})

// for(let i=0; i < 10; i++) console.log(seqResidues[0][i][0].residue)

const chainSequences = seqChains.map((seq, index) => {
    let n = 0
    return (
        seq.map((residue, i) => {            
           if(n <= seqResidues[index].length-1) { 
                if(residue === seqResidues[index][n][0].residue && i === seqResidues[index][n][0].residue_index-1) {
                    n++
                    return seqResidues[index][n-1]
                } else {
                    return null
                }
            } else return null
        })
    )
})

const residueBackbones = chainSequences.map((chain, index) => chain.map((residue, i) => {
    if(!residue) {
        return null
    } else return residue.filter((atom) => atom.atom_type === 'N' || atom.atom_type === 'C' || atom.atom_type === 'CA')
}))
    


// if(index === 0 && i === 10) console.log(residue.filter((atom) => atom.atom_type === 'N' || atom.atom_type === 'C' || atom.atom_type === 'CA'))
    // return
// console.log(partialBackbones[0])
// // residue.filter((atom) => atom.atom_type === 'N' || atom.atom_type === 'CA' || atom.atom_type === 'C')))
// console.log(partialBackbones[0][0])
// sequence with holes end

// view
// sequences[0].forEach((residue, index) => {
//     // if(index === 0 ) console.log(sequences[0][index])
//     if(chainSequences[0][index]) {
//         if(chainSequences[0][index] !== null) {
//             if(residue === chainSequences[0][index][0].residue) {
//                 console.log(residue)
//             }
//         }
//     } else if(chainSequences[0][index] === null) {
//         console.log(null)
//     } else console.log(index)
// })
// view




let chainObject = {}
chainAtoms.forEach((chain, index) => {
    chainObject = {...chainObject, [chainLabels[index]]: chain } })

object = {...object, chains: chainObject}
object = {...object, atoms: newAtoms}


const backbones = chainAtoms.map((chain) => chain.filter((residue) => residue.atom_type === 'CA' || residue.atom_type === 'C' || residue.atom_type === 'N' ))

object = {...object, backbones: backbones }


// torsion angles

const residues = backbones.map((chain) => {
    return (
        chain.filter((atom, index, array) => {
            if(index !== array.length-1) {
            return (
                atom.residue !== array[index+1].residue
            )
        }
    }).map((atom) => atom.residue)
    )
})

const residueAtoms = residues.map((chain, index) => chain.map((residue, i) => backbones[index].slice(i*3, i*3+3)))
const residueAtomsWLabels = residues.map((chain, index) => chain.map((residue, i) => {const resObj = { [residue]: backbones[index].slice(i*3, i*3+3) }; return resObj}))



const torsionAngles = residueBackbones.map((chain, i) => {
    return (
        chain.map((residue, index, array) => {
            // if(index !== 0 && index !== array.length-1) {
                if(residue && array[index-1] && array[index+1]) {
                    const tag = Object.keys(residue)[0]
                    const tagix = Object.keys(array[index-1])[0]
                    const tagii = Object.keys(array[index+1])[0]
                    const vectors = { 
                        Ni: createVectorObj([residue[0].x, residue[0].y, residue[0].z]),

                        Cix: createVectorObj([array[index-1][2].x, array[index-1][2].y, array[index-1][2].z]), 
                        Cia: createVectorObj([residue[1].x, residue[1].y, residue[1].z]),

                        Nii: createVectorObj([array[index+1][0].x, array[index+1][0].y, array[index+1][0].z]), 
                        Ci: createVectorObj([residue[2].x, residue[2].y, residue[2].z]) 
                }


                const phiPlanes = [[vectors.Cix, vectors.Ni, vectors.Cia], [vectors.Ni, vectors.Cia, vectors.Ci]]
                const psiPlanes = [[vectors.Ni, vectors.Cia, vectors.Ci], [vectors.Cia, vectors.Ci, vectors.Nii]]
                
                const phiNormals = phiPlanes.map((plane) => {
                    const U = v.subVector(plane[0], plane[1])
                    const W = v.subVector(plane[0], plane[2])
                    const V = v.crossProduct(U, W)
                    const Vi = v.unitVector(V)

                    return Vi
                })

                const psiNormals = psiPlanes.map((plane) => {
                    const U = v.subVector(plane[0], plane[1])
                    const W = v.subVector(plane[0], plane[2])
                    const V = v.crossProduct(U, W)
                    const Vi = v.unitVector(V)
                    
                    return Vi
                })

                const phiCross = v.crossProduct(phiNormals[0], phiNormals[1])
                const psiCross = v.crossProduct(psiNormals[0], psiNormals[1])

                // const phiSign = Math.sign((phiNormals[1].i - phiNormals[0].i)*(phiNormals[1].j - phiNormals[0].j)*(phiNormals[1].k - phiNormals[0].k))
                // const psiSign = Math.sign((psiNormals[1].i - psiNormals[0].i)*(psiNormals[1].j - psiNormals[0].j)*(psiNormals[1].k - psiNormals[0].k))

                const phiSign = Math.sign(phiCross.k)
                const psiSign = Math.sign(psiCross.k)

                const phi = Math.acos(v.dotProduct(phiNormals[0], phiNormals[1])) * 180/Math.PI*phiSign
                const psi = Math.acos(v.dotProduct(psiNormals[0], psiNormals[1])) * 180/Math.PI*psiSign

                const angles = {phi: phi, psi: psi}

                return angles

            } else return null
        })
    )
})

console.log(residueBackbones[0][2][0].residue)
console.log(torsionAngles[0][2])
console.log(residueBackbones[0][10][0].residue)
console.log(torsionAngles[0][10])
console.log(residueBackbones[0][11][0].residue)
console.log(torsionAngles[0][11])
console.log(residueBackbones[0][59][0].residue)
console.log(torsionAngles[0][59])
console.log(residueBackbones[0][94][0].residue)
console.log(torsionAngles[0][94])

// console.log(torsionAngles[0])

// if( index === 1 && i === 0) {
                    //     const angles = [45, 135, 225, 315]
                    //     const cos = []
                    //     const sin = []
                    //     angles.forEach((angle) => cos.push(Math.cos(angle*Math.PI/180)))
                    //     angles.forEach((angle) => sin.push(Math.sin(angle*Math.PI/180)))
                    //     const vectors = []
                    //     angles.forEach((angle, index) => vectors.push(v.createVectorObj([cos[index], sin[index], 0])))
                    //     const cross = []
                    //     const base = v.createVectorObj([1, 0, 0])
                    //     vectors.forEach((vector) => cross.push(v.crossProduct(vector, base)))
                    //     console.log(vectors)
                    //     console.log(cross)
                    // }

 // if( index === 100 && i === 0) {
                //     console.log(psiNormals)
                //     const sign = (psiNormals[1].i - psiNormals[0].i)*(psiNormals[1].j - psiNormals[0].j)*(psiNormals[1].k - psiNormals[0].k)
                //     console.log(sign)
                //     console.log(Math.sign(sign))
                //     console.log(psi)
                // }
// console.log(residueAtoms[0])

let torsionObject = {}
residues.forEach((chain, index) => { torsionObject = {...torsionObject, chain: chainLabels[index], [chainLabels[index]]: chain.map((residue, i, array) => { if(i !== 0 && i !== array.length-1) {const obj = { [residue]: torsionAngles[index][i]}; return obj } else { const blank = { [residue]: {phi: null, psi: null} }; return blank } })} } )

object = { ...object, torsion_angles: torsionObject }


    // console.log(object.backbones[0].length, object.backbones[1].length, object.backbones[2].length, object.backbones[3].length, object.backbones[4].length, object.backbones[5].length)
    // console.log(object.chains[0].length, object.chains[1].length, object.chains[2].length, object.chains[3].length, object.chains[4].length, object.chains[5].length)


