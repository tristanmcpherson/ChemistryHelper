// Learn more about F# at http://fsharp.org

open System
open FParsec
open FSharp.Data

type Molecule = { Identifier : string; AtomCount : int }

type Compound = { Parts : Molecule list }

type Reactant = { Coefficient : int; Parts : Compound list}

type Reaction = { LHS : Reactant list; RHS : Reactant list }

type PeriodicTable = CsvProvider<"PeriodicTable.csv">
let periodicTable = 
    PeriodicTable.Load("PeriodicTable.csv").Rows |>
    Seq.map (fun r -> (r.Symbol.Trim(), float r.Atomic_Weight)) |>
    Map.ofSeq

let elementNames =
    periodicTable |> Map.toArray |> Array.map (fun (k, _) -> k) |> Set.ofArray

let stringElement s : Parser<string,'u> =
    fun stream ->
        let curr = stream.Read()
        let peek = stream.Peek()
        let first = new String([|curr|])
        let str = new String([|curr; peek|])
        if Set.contains str s then stream.Read() |> ignore; Reply(str)
        else if Set.contains first s then Reply(first) 
             else Reply(Error, NoErrorMessages)

let str = pstring

//let element = manyMinMaxSatisfy2 1 2 Char.IsUpper Char.IsLower
let element = stringElement elementNames
let molecule = (element .>>. opt (str "_" >>. pint32)) |>> (fun (e, c) -> { Identifier = e; AtomCount = Option.defaultValue 0 c})
let compound = (str "(" |> optional) >>. (many molecule) .>> (str ")" |> optional) |>> (fun l -> { Parts = l })
let reactant = opt pint32 .>>. many compound |>> (fun (n, l) -> { Coefficient = Option.defaultValue 1 n; Parts = l})
let reactionPart = sepEndBy1 reactant (str "+")
let reaction = reactionPart .>> (str "->") .>>. reactionPart |>> (fun (l, r) -> { LHS = l; RHS = r })

let printMolecule (m:Molecule) =
    let count = if m.AtomCount = 1 then "" else "_" + string m.AtomCount
    String.Format ("{0}{1}", m.Identifier, count)

let printCompound (c:Compound) =
    List.map printMolecule c.Parts |>
    String.concat ""

let printReactant (r:Reactant) =
    let molecules = String.concat " " (List.map (printCompound) r.Parts)
    let a = if r.Coefficient = 1 then "" else string r.Coefficient
    String.Format("{0}{1}", a, molecules)

let printReaction (r:Reaction) = 
    let reactant = List.map printReactant >>
                   String.concat "+" 
    String.Format ("{0} -> {1}", reactant r.LHS, reactant r.RHS)

[<EntryPoint>]
let main argv =
    let m = match run element "" with
            | Success(result,_,_) -> result
            | Failure(msg,_,_) -> failwith msg
    printfn "%s" <| printReaction m
    //printfn "%s %i" m.Identifier m.AtomCount
    printfn "Hello World from F#!"
    0 // return an integer exit code
