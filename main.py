import Bio.PDB
import numpy as np
import matplotlib.pyplot as plt


def find_angle(p1, p2, p3, p4):
	q1 = p2 - p1
	q2 = p3 - p2
	q3 = p4 - p3

	n1 = np.cross(q1, q2)
	n2 = np.cross(q2, q3)

	n1 = n1 / np.sqrt(np.sum(n1 * n1))
	n2 = n2 / np.sqrt(np.sum(n2 * n2))

	q2 = q2 / np.sqrt(np.sum(q2 * q2))
	n3 = np.cross(n1, n2)
	n3 = n3 / np.sqrt(np.sum(n3 * n3))

	angle = np.arccos(np.sum(n1 * n2)) * np.sign(np.sum(q2 * n3))

	return angle


if __name__ == "__main__":
	PDB_NAMES = ('4d2i', '4b2i', '3a2b', '4q4w')

	for pdb_file in PDB_NAMES:
		for model in Bio.PDB.PDBParser().get_structure(pdb_file, pdb_file + '.pdb'):
			for chain in model:
				phi_angles = []
				psi_angles = []

				name = chain.get_full_id()
				name_str = ''
				for clever_stuff in name:
					name_str += str(clever_stuff) + '_'

				poly = Bio.PDB.Polypeptide.Polypeptide(chain)

				poly_len = len(chain)
				prev_residue = None

				for i, residue in enumerate(chain):
					try:
						N = residue['N'].coord
						CA = residue['CA'].coord
						C = residue['C'].coord
					# If no element.
					except Exception:
						continue

					if i > 0:
						try:
							prev_C = prev_residue['C'].coord
						except Exception:
							continue

						phi = find_angle(prev_C, N, CA, C)
						phi_angles.append(phi)
					else:
						phi_angles.append(None)

					if i > 0:
						try:
							prev_N = prev_residue['N'].coord
							prev_CA = prev_residue['CA'].coord
							prev_C = prev_residue['C'].coord
						except Exception:
							continue

						psi = find_angle(prev_N, prev_CA, prev_C, N)
						psi_angles.append(psi)

					prev_residue = residue

				psi_angles.append(None)

				my_angles = np.stack((phi_angles, psi_angles), axis=1)

				plt.scatter(phi_angles, psi_angles)
				plt.savefig(name_str + '.png')
				plt.close()
