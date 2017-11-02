# Transcranial direct current stimulation over the left prefrontal cortex increases randomness of choice in instrumental learning

todo: [![DOI](https://zenodo.org/badge/19634/ihrke/2016-placebo-tdcs-study.svg)](https://zenodo.org/badge/latestdoi/19634/ihrke/2016-placebo-tdcs-study) 

This repository contains data for the paper "Transcranial direct current stimulation over the left prefrontal cortex increases randomness of choice in instrumental learning".

If you want to use this data in a research publication,
please cite [our paper](https://www.sciencedirect.com/science/article/pii/S0010945214002895).


Turi, Z., Mittner, M., Opitz, A., Popkes, M., Paulus, W. & Antal, A. (2015).
Transcranial direct current stimulation over the left prefrontal cortex increases randomness of choice in instrumental learning. Cortex. 63:145-154. doi:10.1016/j.cortex.2014.08.026

~~~{bibtex}
@article{Turi_tdcs2015,
  author={Turi, Z. and Mittner, M. and Opitz, A. and Popkes, M. and Paulus, W. and Antal, A.},
  title={Transcranial direct current stimulation over the left prefrontal cortex increases randomness of choice in instrumental learning},
  year=2015,
  journal={Cortex},
  volume=63,
  number=145-154,
  doi=10.1016/j.cortex.2014.08.026
}
~~~

## Structure of the paradigm in brief:

- Learning phase: AB, CD and EF symbol pairs presented 120 times each, 360 times in total
- Transfer phase: all the possible item combinations (e.g., AB, AC, AE, AF...) presented 12 times

## Stimulation code:
- A = real tDCS
- B = sham tDCS

## Data

- Raw data for the learning phase is located in `data/raw/Learning`
- Raw data for the transfer phase is located in `data/raw/Transfer`
- Subfolders A and B code the type of stimulation 
- Raw data is provided in `.txt` format

## Output for the LEARNING phase:

Variables are coded as follows:

- `Block` - block number (1, 2, 3, 4, 5, 6), each of the 3 symbol pairs were presented 20 times/block
- `Pair_type`  - symbol pairs with 80/20% (AB), 70/30% (CD) or 60/40% (EF) reward contingency
- `Feedback_type` - type of feedback: 1 reward was received, 0 no reward
- `Side` - presentation of the items on the left (1) and right (2) side
- `ACC` - accuracy: 1 items with the higher reward probability were chosen (e.g., A over B), 0 items with lower rew.prob. were chosen (e.g., B over A), missed no decision was made within 1700 msec
- `RT` - reaction time in msec

## Output for the TRANSFER phase:

Variables are coded as follows:

- `Side` - presentation of the items on the left (1) and right (2) side
- `RT` - reaction time in msec
- `ACC` - accuracy: 1 items with the higher reward probability were chosen (e.g., A over B), 0 items with lower rew.prob. were chosen (e.g., B over A), missed no decision was made within 1700 msec
- `Type`  - code for all possible item combinations 
	- 1:AB
	- 2:AC
	- 3:AD
	- 4:AE
	- 5:AF
	- 6:BC
	- 7:BD
	- 8:BE
	- 9:BF
	- 10:CD
	- 11:CE
	- 12:CF
	- 13:DE
	- 14:DF
	- 15:EF