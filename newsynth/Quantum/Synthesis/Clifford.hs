{-# OPTIONS_GHC -fno-warn-incomplete-patterns #-}

-- | This module provides an efficient symbolic representation of the
-- Clifford group on one qubit. This group is generated by /S/, /H/,
-- and the scalar ω = [exp /i/π\/4]. It has 192 elements. 

module Quantum.Synthesis.Clifford (
  -- * The Clifford group
  Clifford,
  
  -- ** Constructors
  clifford_X,
  clifford_Y,
  clifford_Z,
  clifford_H,
  clifford_S,
  clifford_SH,
  clifford_E,
  clifford_W,
  ToClifford(to_clifford),
  
  -- ** Deconstructors
  clifford_decompose,
  Axis(..),
  clifford_decompose_coset,
  
  -- ** Group operations
  clifford_id,
  clifford_mult,
  clifford_inv,
  
  -- ** Conjugation by /T/
  clifford_tconj  
  ) where

-- ----------------------------------------------------------------------
-- * The Clifford group

-- $ We could, in principle, implement the Clifford group as an
-- enumerated type with 192 elements, and a large 192×192 lookup
-- table for the group multiplication. Instead, we take advantage of
-- some of the internal structure of the group to reduce the size of
-- the lookup tables. The resulting implementation is still very
-- efficient.

-- | A type representing single-qubit Clifford operators.
data Clifford = Clifford Int Int Int Int
                deriving (Eq, Ord)

instance Show Clifford where
  show (Clifford a b c d) = "C" ++ show a ++ show b ++ show c ++ show d

-- ----------------------------------------------------------------------
-- ** Constructors

-- | The Pauli /X/-gate as a Clifford operator.
clifford_X :: Clifford
clifford_X = Clifford 0 1 0 0

-- | The Pauli /Y/-gate as a Clifford operator.
clifford_Y :: Clifford
clifford_Y = Clifford 0 1 2 2

-- | The Pauli /Z/-gate as a Clifford operator.
clifford_Z :: Clifford
clifford_Z = Clifford 0 0 2 0

-- | The Hadamard gate as a Clifford operator.
clifford_H :: Clifford
clifford_H = Clifford 1 0 1 5

-- | The Clifford operator /S/.
clifford_S :: Clifford
clifford_S = Clifford 0 0 1 0

-- | The Clifford operator /SH/.
clifford_SH :: Clifford
clifford_SH = clifford_S `clifford_mult` clifford_H

-- | The Clifford operator /E/ = /H//S/[sup 3]ω[sup 3]. This operator is
-- uniquely determined by the properties /E/³ = /I/, 
-- /EXE/⁻¹ = /Y/, /EYE/⁻¹ = /Z/, and /EZE/⁻¹ = /X/.
-- 
-- \[image E.png]
clifford_E :: Clifford
clifford_E = Clifford 1 0 0 0

-- | The Clifford operator ω = [exp /i/π\/4].
clifford_W :: Clifford
clifford_W = Clifford 0 0 0 1

-- | A type class for things that can be exactly converted to a
-- Clifford operator. One particular instance of this is 'String', so
-- that Clifford operators can be denoted, e.g.,
-- 
-- > to_clifford "-iX"
-- 
-- The valid characters for such string conversions are @\"XYZHSEIWi-\"@.
class ToClifford a where
  -- | Convert any suitable thing to a Clifford operator.
  to_clifford :: a -> Clifford
  
instance ToClifford Clifford where
  to_clifford = id
  
instance ToClifford Char where
  to_clifford 'E' = clifford_E
  to_clifford 'X' = clifford_X
  to_clifford 'S' = clifford_S
  to_clifford 'W' = clifford_W
  to_clifford 'I' = clifford_id
  to_clifford 'i' = Clifford 0 0 0 2
  to_clifford '-' = Clifford 0 0 0 4
  to_clifford 'H' = clifford_H
  to_clifford 'Y' = clifford_Y
  to_clifford 'Z' = clifford_Z
  to_clifford x = error $ "ToClifford Char: unknown gate " ++ show x

instance ToClifford a => ToClifford [a] where
  to_clifford [] = clifford_id
  to_clifford (h:t) = to_clifford h `clifford_mult` to_clifford t

-- ----------------------------------------------------------------------
-- ** Deconstructors

-- | Given a Clifford operator /U/, return (/a/, /b/, /c/, /d/) such that
-- 
-- * /U/ = /E/[sup /a/]/X/[sup /b/]/S/[sup /c/]ω[sup /d/],
-- 
-- * /a/ ∈ {0, 1, 2}, /b/ ∈ {0, 1}, /c/ ∈ {0, …, 3}, and /d/ ∈ {0, …,
-- 7}.
-- 
-- Here, /E/ = /H//S/[sup 3]ω[sup 3]. Note that /E/, /X/, /S/, and ω have order
-- 3, 2, 4, and 8, respectively. Moreover, each Clifford operator can
-- be uniquely represented as above.
clifford_decompose :: (ToClifford a) => a -> (Int, Int, Int, Int)
clifford_decompose m = (a,b,c,d) where
  Clifford a b c d = to_clifford m

-- | A axis is either /I/, /H/, or /SH/.
data Axis = Axis_I | Axis_H | Axis_SH
           deriving (Eq, Show)

instance ToClifford Axis where
  to_clifford Axis_I = to_clifford "I"
  to_clifford Axis_H = to_clifford "H"
  to_clifford Axis_SH = to_clifford "SH"

-- | Given a Clifford operator /U/, return (/K/, /b/, /c/, /d/) such that
-- 
-- * /U/ = /K//X/[sup /b/]/S/[sup /c/]ω[sup /d/],
-- 
-- * /K/ ∈ {/I/, /H/, /SH/}, /b/ ∈ {0, 1}, /c/ ∈ {0, …, 3}, and /d/ ∈ {0, …,
-- 7}.
clifford_decompose_coset :: (ToClifford a) => a -> (Axis, Int, Int, Int)
clifford_decompose_coset u = case op of
  Clifford 0 b c d -> (Axis_I, b, c, d)
  Clifford 1 b c d -> (Axis_H, b', c', d') where
    Clifford 0 b' c' d' = clifford_inv "H" `clifford_mult` op
  Clifford 2 b c d -> (Axis_SH, b', c', d') where
    Clifford 0 b' c' d' = clifford_inv "SH" `clifford_mult` op
  where
    op = to_clifford u
  
-- ----------------------------------------------------------------------
-- ** Group operations

-- | The identity Clifford operator.
clifford_id :: Clifford
clifford_id = Clifford 0 0 0 0

-- | Clifford multiplication.
clifford_mult :: Clifford -> Clifford -> Clifford
clifford_mult u1 u2 = u where
  -- U = U1 U2   
  --   = A1 B1 C1 D1 A2 B2 C2 D2
  --   = A1 (B1 C1 A2) B2 C2 D1 D2
  --   = A1 (A3 B3 C3 D3) B2 C2 D1 D2
  --   = A1 A3 B3 (C3 B2) C2 D3 D1 D2
  --   = A1 A3 B3 (B2 C4 D4) C2 D3 D1 D2
  --   = (A1 A3) (B3 B2) (C4 C2) (D4 D3 D1 D2)
  --   = A B C D
  Clifford a1 b1 c1 d1 = u1
  Clifford a2 b2 c2 d2 = u2
  (a3, b3, c3, d3) = conj3 b1 c1 a2
  (c4, d4) = conj2 c3 b2
  a = (a1 + a3) `mod` 3
  b = (b3 + b2) `mod` 2
  c = (c4 + c2) `mod` 4
  d = (d4 + d3 + d1 + d2) `mod` 8
  u = Clifford a b c d

-- | Clifford inverse.
clifford_inv :: (ToClifford a) => a -> Clifford
clifford_inv op = Clifford a2 b2 c2 d3 where
  -- U⁻¹ = (A B C)⁻¹ D⁻¹ = (A2 B2 C2 D2) D⁻¹
  Clifford a b c d = to_clifford op
  (a2, b2, c2, d2) = cinv a b c
  d3 = (d2 - d) `mod` 8

-- ----------------------------------------------------------------------
-- ** Conjugation by /T/

-- | Given a Clifford gate /C/, return an axis /K/ ∈ {/I/, /H/, /SH/}
-- and a Clifford gate /C'/ such that
-- 
-- * /C//T/ = /K//T//C/'.
clifford_tconj ::  Clifford -> (Axis, Clifford)
clifford_tconj u = (k, v) where
  -- U T = A1 B1 C1 D1 T
  --     = (A1 B1 T) C1 D1
  --     = (K T B1 C2 D2) C1 D1
  --     = K T B1 (C2 C1) (D2 D1)
  Clifford a1 b1 c1 d1 = u
  (k, c2, d2) = tconj a1 b1
  c = (c2 + c1) `mod` 4
  d = (d2 + d1) `mod` 8
  v = Clifford 0 b1 c d

-- ----------------------------------------------------------------------
-- ** Lookup tables

-- | 'conj2' /c/ /b/ returns (/c/', /d/') such that
-- 
-- * /S/[sup /c/]/X/[sup /b/] = /X/[sup /b/]/S/[sup /c/']ω[sup /d/'].
conj2 :: Int -> Int -> (Int, Int)
conj2 0 0 = (0,0)
conj2 0 1 = (0,0)
conj2 1 0 = (1,0)
conj2 1 1 = (3,2)
conj2 2 0 = (2,0)
conj2 2 1 = (2,4)
conj2 3 0 = (3,0)
conj2 3 1 = (1,6)

-- | 'conj3' /b/ /c/ /a/ returns (/a/', /b/', /c/', /d/') such that
-- 
-- * /X/[sup /b/]/S/[sup /c/]/E/[sup /a/] = /E/[sup /a/']/X/[sup /b/']/S/[sup /c/']ω[sup /d/'].
conj3 :: Int -> Int -> Int -> (Int, Int, Int, Int)
conj3 0 0 0 = (0,0,0,0)
conj3 0 0 1 = (1,0,0,0)
conj3 0 0 2 = (2,0,0,0)
conj3 0 1 0 = (0,0,1,0)
conj3 0 1 1 = (2,0,3,6)
conj3 0 1 2 = (1,1,3,4)
conj3 0 2 0 = (0,0,2,0)
conj3 0 2 1 = (1,1,2,2)
conj3 0 2 2 = (2,1,0,0)
conj3 0 3 0 = (0,0,3,0)
conj3 0 3 1 = (2,1,3,6)
conj3 0 3 2 = (1,0,1,2)
conj3 1 0 0 = (0,1,0,0)
conj3 1 0 1 = (1,0,2,0)
conj3 1 0 2 = (2,1,2,2)
conj3 1 1 0 = (0,1,1,0)
conj3 1 1 1 = (2,1,1,0)
conj3 1 1 2 = (1,1,1,0)
conj3 1 2 0 = (0,1,2,0)
conj3 1 2 1 = (1,1,0,6)
conj3 1 2 2 = (2,0,2,6)
conj3 1 3 0 = (0,1,3,0)
conj3 1 3 1 = (2,0,1,4)
conj3 1 3 2 = (1,0,3,2)

-- | 'cinv' /a/ /b/ /c/ returns (/a/', /b/', /c/', /d/') such that
-- 
-- * (/E/[sup /a/]/X/[sup /b/]/S/[sup /c/])⁻¹ = /E/[sup /a/']/X/[sup /b/']/S/[sup /c/']ω[sup /d/'].
cinv :: Int -> Int -> Int -> (Int, Int, Int, Int)
cinv 0 0 0 = (0,0,0,0)
cinv 0 0 1 = (0,0,3,0)
cinv 0 0 2 = (0,0,2,0)
cinv 0 0 3 = (0,0,1,0)
cinv 0 1 0 = (0,1,0,0)
cinv 0 1 1 = (0,1,1,6)
cinv 0 1 2 = (0,1,2,4)
cinv 0 1 3 = (0,1,3,2)
cinv 1 0 0 = (2,0,0,0)
cinv 1 0 1 = (1,0,1,2)
cinv 1 0 2 = (2,1,0,0)
cinv 1 0 3 = (1,1,3,4)
cinv 1 1 0 = (2,1,2,2)
cinv 1 1 1 = (1,1,1,6)
cinv 1 1 2 = (2,0,2,2)
cinv 1 1 3 = (1,0,3,4)
cinv 2 0 0 = (1,0,0,0)
cinv 2 0 1 = (2,1,3,6)
cinv 2 0 2 = (1,1,2,2)
cinv 2 0 3 = (2,0,3,6)
cinv 2 1 0 = (1,0,2,0)
cinv 2 1 1 = (2,1,1,6)
cinv 2 1 2 = (1,1,0,2)
cinv 2 1 3 = (2,0,1,6)

-- | 'tconj' /a/ /b/ returns (/K/, /c/, /d/) such that
-- 
-- * /E/[sup /a/]/X/[sup /b/]/T/ = /K//T//X/[sup /b/]/S/[sup /c/]ω[sup /d/].
tconj :: Int -> Int -> (Axis, Int, Int)
tconj 0 0 = (Axis_I,  0, 0)
tconj 0 1 = (Axis_I,  1, 7)
tconj 1 0 = (Axis_H,  3, 3)
tconj 1 1 = (Axis_H,  2, 0)
tconj 2 0 = (Axis_SH, 0, 5)
tconj 2 1 = (Axis_SH, 1, 4)
