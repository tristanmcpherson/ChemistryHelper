// Learn more about F# at http://fsharp.org

open System
open Microsoft.FSharp.Core
open FSharp.Data
open FsAlg.Generic

module List =
    let public pairwiseFold (folder : 'State -> 'T -> 'T -> 'State) (state : 'State) (source : 'T list) =
        List.pairwise source |>
        List.fold (fun a (f, s) -> folder a f s) state
           
    let public findIndices (pred : 'a -> bool) (chars : 'a list) =
        List.indexed chars |>
        List.filter (fun (_, e) -> pred e) |>
        List.map (fun (i, _) -> i)

    let public toTuple (a : 'a list) =
        match a with
        | [a; b] -> Some(a, b)
        | _ -> None

module String =
    let public toChars (str : string) =
        str.ToCharArray() |> List.ofArray

    let public splitAt (str : string) (delim : string) = 
        List.ofArray (str.Split(delim)) |> 
        List.toTuple

    let public splitAtBack (str : string) (delim : char) =
        let chars = toChars str
        List.findIndexBack ((=) delim) chars |>
        (fun i -> List.splitAt i chars) |>
        (fun (f, s) -> (f, List.skip 1 s))

    let public splitBy (pred : char -> bool) (s : string) =
        List.findIndices pred (toChars s) @ [s.Length] |>
        List.pairwiseFold (fun a e e1 -> s.Substring(e, e1 - e) :: a) (List.empty) |>
        List.rev

type PeriodicTable = CsvProvider<"PeriodicTable.csv">
let periodicTable = 
    let csv = PeriodicTable.Load("PeriodicTable.csv")
    csv.Rows |>
    Seq.map (fun r -> (r.Symbol.Trim(), float r.Atomic_Weight)) |>
    Map.ofSeq

type Molecule = { Identifier : string; AtomCount : int }

type Reactant = { Coefficient : int; Parts : Molecule list}

type Reaction = { LHS : Reactant list; RHS : Reactant list }

let splitString s d =
    match String.splitAt s d with
    | Some (a, b) -> (a, b)
    | None -> failwithf "Failed to split string %s at %s" s d

let printMolecule (m:Molecule) =
    let count = if m.AtomCount = 1 then "" else "_" + string m.AtomCount
    String.Format ("{0}{1}", m.Identifier, count)

let printReactant (r:Reactant) =
    let molecules = String.concat " " (List.map (printMolecule) r.Parts)
    let a = if r.Coefficient = 1 then "" else string r.Coefficient
    String.Format("{0}{1}", a, molecules)

let printReaction (r:Reaction) = 
    let reactant = List.map printReactant >>
                   String.concat "+" 
    String.Format ("{0} -> {1}", reactant r.LHS, reactant r.RHS)

let parseMolecule (s:string) : Molecule = 
    let split = String.splitAt s "_"
    match split with 

    | Some (m, c) -> { Identifier = m; AtomCount = int c }
    | None -> { Identifier = s; AtomCount = 1 }

let parseReactant (s:string) : Reactant =
    let chars = s |> String.toChars
    let stringify (c : char list) =
        c |> Array.ofList |> String.Concat

    let uppers = chars |> List.findIndices (Char.IsUpper)
    let lp = chars |> List.findIndices ((=) '(')
    let rp = chars |> List.findIndices ((=) ')') |> List.map (fun x -> (uppers |> List.find (fun y -> y > x)))

    let parse s =
        List.map parseMolecule (String.splitBy (fun x -> Char.IsUpper x) s)

    let groups = List.zip lp rp |> 
                 List.map (fun (st, e) -> s.Substring(st, e - st)) |>
                 List.map (fun (x : string) -> x.Replace(")", "").Replace("(", "")) |>
                 List.map (String.toChars) |>
                 List.map (fun x -> String.splitAtBack (stringify x) '_') |>
                 List.map (fun (d, b) -> stringify (List.concat [['(']; d; [')']; b]), parse (stringify d), stringify b |> int) |>
                 List.map (fun (o, ms, n) -> (o, List.fold (fun (a:string) (m : Molecule) -> String.Concat [a; m.Identifier; "_"; string (m.AtomCount * n)]) "" ms))

    let resolved = List.fold (fun (a : string) ((o : string), (r : string)) -> a.Replace(o, r)) s groups

    let start = if lp.Length = 0 then uppers.Head else if uppers.Head < lp.Head then uppers.Head else lp.Head
    let split = List.splitAt start (String.toChars resolved)

    match split with
    | [], a -> { Coefficient = 1;                     Parts = a |> stringify |> parse }
    | c, a  -> { Coefficient = c |> stringify |> int; Parts = a |> stringify |> parse }

let parseReaction (s:string) : Reaction =
    let (l, r) = splitString s "->"
    let getReactants (x:string) = x.Split "+" |>
                                  (fun m -> List.map parseReactant (List.ofArray m))
                         
    { LHS = getReactants l; RHS = getReactants r; }

let balanceEquation (reaction : Reaction) = 
    let elementLookup = List.collect (fun (x : Reactant) -> List.map (fun (y : Molecule) -> y.Identifier) x.Parts) >> 
                        List.distinct
    let finalLookup = (elementLookup reaction.LHS) @ (elementLookup reaction.RHS) |> 
                      List.distinct |> 
                      List.map (fun x -> (x, 0)) |> Map.ofList

    let getVector = List.fold (fun a (x : Molecule) -> Map.add x.Identifier (a.[x.Identifier] + x.AtomCount) a) finalLookup
    let getVectors = List.map (fun (x : Reactant) -> getVector x.Parts)

    let vectorize s = getVectors >>
                      List.map (fun x -> Map.toArray x) >>
                      List.map (fun x -> Array.map (fun (_, b) -> float (s b)) x) >>
                      List.map (fun x -> Vector x)

    let vectors = (vectorize (~+) reaction.LHS) @ (vectorize (~-) reaction.RHS)

    let rev = List.rev vectors
    let balanced = Matrix.ofCols rev.Tail |> (fun m -> Matrix.solve m (rev.Head))
    balanced |> (fun v -> Vector (Array.map (~-) (v.ToArray())))

let convertUnits f t value =
    let map = [("kg", "g"), 1E3;
               ("g", "mg"), 1E3;] |> Map.ofList

    let find k v k' v' = 
        match k, v with 
        | k, v when k = (k', v') -> Some (v, ( * ))
        | k, v when k = (v', k') -> Some (v, ( / ))
        | _ -> None

    match Map.tryPick (fun k v -> find k v f t) map with 
     | Some (v, o) -> o value v
     | None -> failwith "No conversion found"

let molecularMass =
    (fun m -> periodicTable.[m.Identifier] * float m.AtomCount)

let molConversion reactant mass =
    let gPerMole = List.map (molecularMass) reactant.Parts |>
                   List.sum

    (mass, gPerMole, mass / gPerMole / float reactant.Coefficient)

let parseInput (s : string) =
    let pair = splitString s " "
    pair |> (fun (a, b) -> if b <> "g" then (convertUnits b "g" (float a)) else float a)

    //let getVector = List.countBy (fun (x : Molecule) -> x.Identifier) >> Map.ofList
    //let getVectors = List.map (fun (x : Reactant) -> getVector x.Parts)
    //getVectors reaction.LHS |> List.map FsAlg.Generic.Vector
    
[<EntryPoint>]
let main argv =
    let input = "C_4H_8+O_2->CO_2+H_2O"
    let reaction = parseReaction input
    let mass = convertUnits "kg" "g" 5.22
    printfn "Parsed equation as: %s" <| printReaction reaction
    let balanced = balanceEquation reaction
    printfn "Mass of reactant: %fg, molar mass: %fg/mol, mol: %f" <||| molConversion reaction.LHS.Head mass
    0 // return an integer exit code
