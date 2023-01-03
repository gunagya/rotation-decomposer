{-# LANGUAGE ForeignFunctionInterface #-}

module LsqeccGlue.ExportCCalls where

import System.IO
import System.Random

import Foreign.C.Types
import Foreign.C.String
import Foreign.Ptr
import Foreign.Storable

import Quantum.Synthesis.GridSynth
import Quantum.Synthesis.CliffordT
import Quantum.Synthesis.SymReal

k_defaultEffort = 25 -- Value recommended in the gridsynth source

nullChar :: CChar
nullChar = fromIntegral 0


type CIntExitCode = CInt
k_exitCodeSuccess :: CIntExitCode
k_exitCodeSuccess = 0

k_exitCodeBufferTooSmall :: CIntExitCode
k_exitCodeBufferTooSmall = 2

k_angleParseError :: CIntExitCode
k_angleParseError = 3


writeGates :: String -> CInt -> Ptr CChar -> IO CInt
writeGates _ (0) _ = return k_exitCodeBufferTooSmall
writeGates [] _ ptr = do
    poke ptr nullChar
    return k_exitCodeSuccess
writeGates (g:gs) n ptr = do
    poke ptr (castCharToCChar g) -- Write g into the adress at ptr
    writeGates gs (n-1) (plusPtr ptr 1)


-- Exit codes: 0 = success. 1 = Buffer too small
gridsynth_angle_to_seq_with_precision_hs :: CDouble -> CString -> Ptr CChar -> CInt -> IO CIntExitCode
gridsynth_angle_to_seq_with_precision_hs precision angle_as_cstring buffer buffer_max_bytes = do
    let randomGen = mkStdGen 30 -- The global std random generators induce a crash, so we use the 
    hs_string_angle <- peekCAString angle_as_cstring
    return 0
    case parse_SymReal hs_string_angle of
        Just angle -> do
            let gates = gridsynth_gates randomGen (fromRational $ toRational precision) angle k_defaultEffort
            writeGates (from_gates gates) buffer_max_bytes buffer
        Nothing -> return k_angleParseError
 

 -- TODO: 
 --  1. replace the recursive write gates with one based on forM and offsets
 --  2. Expose effort too?
     



foreign export ccall gridsynth_angle_to_seq_with_precision_hs :: CDouble -> CString -> Ptr CChar -> CInt -> IO CIntExitCode
