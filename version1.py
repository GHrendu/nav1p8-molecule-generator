"""
Nav1.7 Conformation-Aware Virtual Screening
Version 1.0

功能：

1. 读取多个 Nav1.7 构象
2. 读取小分子库
3. 对小分子进行基础过滤
4. 将小分子转换成 docking 所需格式
5. 对每一个 Nav1.7 构象进行 docking
6. 输出 compound × conformation binding matrix
7. 为后续 AI 模型准备数据

注意：
第一版不训练大型 AI 模型。
先把 docking 数据流水线跑通。
"""

import os
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


# ============================================================
# 1. 项目参数
# ============================================================

# Nav1.7 蛋白构象
PROTEIN_STATES = {

    "state_1": "data/proteins/state_1.pdb",
    "state_2": "data/proteins/state_2.pdb",
    "state_3": "data/proteins/state_3.pdb",
    "state_4": "data/proteins/state_4.pdb",

}

# 小分子库
LIGAND_LIBRARY = "data/ligands/library.sdf"

# docking 输出目录
OUTPUT_DIR = "data/results"

# AutoDock Vina 程序位置
########## TODO ##########
VINA_PATH = "vina"


# ============================================================
# 2. docking box 参数
# ============================================================

"""
这里非常重要。

x_center / y_center / z_center
决定 docking 在 Nav1.7 哪个区域进行。

x_size / y_size / z_size
决定 docking box 大小。

这些数值需要根据你的 Nav1.7 binding pocket 来确定。
"""

DOCKING_BOX = {

    "center_x": 0.0,
    "center_y": 0.0,
    "center_z": 0.0,

    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0,

}

########## TODO ##########
# 上面的 docking box 必须根据你最终选择的 Nav1.7 pocket 修改。


# ============================================================
# 3. 小分子过滤参数
# ============================================================

MAX_MW = 500
MAX_LOGP = 5
MAX_TPSA = 140

MIN_MW = 150


def filter_molecule(mol):

    """
    对小分子进行基础 drug-like filtering。

    返回：
        True  = 保留
        False = 删除
    """

    if mol is None:
        return False

    mw = Descriptors.MolWt(mol)

    logp = Crippen.MolLogP(mol)

    tpsa = Descriptors.TPSA(mol)

    if mw < MIN_MW:
        return False

    if mw > MAX_MW:
        return False

    if logp > MAX_LOGP:
        return False

    if tpsa > MAX_TPSA:
        return False

    return True


# ============================================================
# 4. 读取小分子库
# ============================================================

def load_ligands(sdf_file):

    """
    从 SDF 文件读取小分子。

    返回：
        ligand_list

    每个元素：

        {
            "id": molecule ID,
            "mol": RDKit molecule
        }
    """

    supplier = Chem.SDMolSupplier(
        sdf_file,
        removeHs=False
    )

    ligand_list = []

    for i, mol in enumerate(supplier):

        if mol is None:
            continue

        if not filter_molecule(mol):
            continue

        # 尝试读取 molecule ID
        if mol.HasProp("_Name"):
            mol_id = mol.GetProp("_Name")
        else:
            mol_id = f"compound_{i}"

        ligand_list.append(
            {
                "id": mol_id,
                "mol": mol
            }
        )

    print(
        f"原始分子数量：{len(supplier)}"
    )

    print(
        f"过滤后分子数量：{len(ligand_list)}"
    )

    return ligand_list


# ============================================================
# 5. 保存过滤后的 ligand
# ============================================================

def save_filtered_ligands(
        ligand_list,
        output_file
    ):

    """
    将过滤后的分子保存成新的 SDF。
    """

    writer = Chem.SDWriter(output_file)

    for item in ligand_list:

        mol = item["mol"]

        mol.SetProp(
            "compound_id",
            item["id"]
        )

        writer.write(mol)

    writer.close()


# ============================================================
# 6. 运行 AutoDock Vina
# ============================================================

def run_vina(
        receptor,
        ligand,
        output_pose,
        log_file
    ):

    """
    调用 AutoDock Vina。

    receptor:
        Nav1.7 receptor

    ligand:
        ligand pdbqt

    output_pose:
        docking pose 输出文件

    log_file:
        docking log
    """

    command = [

        VINA_PATH,

        "--receptor",
        receptor,

        "--ligand",
        ligand,

        "--center_x",
        str(DOCKING_BOX["center_x"]),

        "--center_y",
        str(DOCKING_BOX["center_y"]),

        "--center_z",
        str(DOCKING_BOX["center_z"]),

        "--size_x",
        str(DOCKING_BOX["size_x"]),

        "--size_y",
        str(DOCKING_BOX["size_y"]),

        "--size_z",
        str(DOCKING_BOX["size_z"]),

        "--out",
        output_pose,

        "--log",
        log_file,

        "--exhaustiveness",
        "8",

        "--num_modes",
        "3",

    ]

    print(
        "Running:",
        " ".join(command)
    )

    subprocess.run(
        command,
        check=True
    )


# ============================================================
# 7. 从 Vina log 中读取 docking score
# ============================================================

def parse_vina_score(log_file):

    """
    从 AutoDock Vina log 文件读取最佳 binding score。

    Vina 输出通常类似：

    mode | affinity | dist from best mode
    1    | -8.4     | ...
    2    | -7.9     | ...

    我们只取第一名。

    返回：

        -8.4
    """

    if not os.path.exists(log_file):
        return None

    with open(
        log_file,
        "r"
    ) as f:

        lines = f.readlines()

    for line in lines:

        line = line.strip()

        if line.startswith("1"):

            parts = line.split()

            try:

                score = float(parts[1])

                return score

            except:

                pass

    return None


# ============================================================
# 8. 建立 docking 数据
# ============================================================

def perform_screening(
        ligand_list
    ):

    """
    对：

        molecule × Nav1.7 state

    进行 docking。

    最终得到：

        [
            {
                compound_id: xxx,
                state: state_1,
                score: -8.4
            },

            {
                compound_id: xxx,
                state: state_2,
                score: -7.8
            }
        ]

    """

    results = []

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True
    )

    for state_name, receptor in PROTEIN_STATES.items():

        print(
            "\n=============================="
        )

        print(
            f"正在筛选 Nav1.7 {state_name}"
        )

        print(
            "=============================="
        )

        state_dir = os.path.join(
            OUTPUT_DIR,
            state_name
        )

        os.makedirs(
            state_dir,
            exist_ok=True
        )

        for item in ligand_list:

            compound_id = item["id"]

            ligand_file = os.path.join(
                state_dir,
                f"{compound_id}.pdbqt"
            )

            pose_file = os.path.join(
                state_dir,
                f"{compound_id}_out.pdbqt"
            )

            log_file = os.path.join(
                state_dir,
                f"{compound_id}.log"
            )

            # ==========================================
            # TODO
            # 这里目前假设 ligand 已经是 pdbqt。
            #
            # 实际使用时，需要增加：
            #
            # SDF
            # ↓
            # 3D structure
            # ↓
            # protonation
            # ↓
            # pdbqt
            #
            # 这一部分可以用 Meeko 完成。
            # ==========================================

            if not os.path.exists(ligand_file):

                print(
                    f"[WARNING] {compound_id} 没有 pdbqt 文件"
                )

                continue

            try:

                run_vina(
                    receptor,
                    ligand_file,
                    pose_file,
                    log_file
                )

                score = parse_vina_score(
                    log_file
                )

            except Exception as e:

                print(
                    f"Docking failed: {compound_id}"
                )

                print(e)

                score = None

            results.append(
                {
                    "compound_id": compound_id,
                    "state": state_name,
                    "score": score
                }
            )

    return pd.DataFrame(results)


# ============================================================
# 9. 生成 Conformational Binding Matrix
# ============================================================

def build_binding_matrix(
        docking_results
    ):

    """
    将：

    compound | state | score

    转换成：

              state_1 state_2 state_3 state_4

    compound1  -8.4   -7.2   -9.1   -8.8

    compound2  -7.1   -9.5   -7.3   -7.4

    """

    matrix = docking_results.pivot(
        index="compound_id",
        columns="state",
        values="score"
    )

    return matrix


# ============================================================
# 10. 计算构象选择性
# ============================================================

def calculate_conformation_selectivity(
        matrix
    ):

    """
    一个非常重要的指标。

    如果：

        State 1 = -9.5
        State 2 = -7.2
        State 3 = -7.0
        State 4 = -7.1

    那么这个分子可能具有：

        State 1 preference

    我们可以计算：

        best_score
        mean_score
        score_range

    """

    output = matrix.copy()

    output["best_score"] = matrix.min(
        axis=1
    )

    output["worst_score"] = matrix.max(
        axis=1
    )

    output["score_range"] = (
        output["worst_score"]
        -
        output["best_score"]
    )

    output["mean_score"] = matrix.mean(
        axis=1
    )

    return output


# ============================================================
# 11. 主程序
# ============================================================

def main():

    print(
        "======================================"
    )

    print(
        "Nav1.7 AI Virtual Screening V1.0"
    )

    print(
        "======================================"
    )

    # --------------------------------------------------------
    # Step 1
    # --------------------------------------------------------

    print(
        "\n[1] Loading ligand library..."
    )

    ligands = load_ligands(
        LIGAND_LIBRARY
    )

    # --------------------------------------------------------
    # Step 2
    # --------------------------------------------------------

    filtered_file = (
        "data/ligands/"
        "filtered_library.sdf"
    )

    save_filtered_ligands(
        ligands,
        filtered_file
    )

    print(
        "\n过滤后的 ligand 已保存：",
        filtered_file
    )

    # --------------------------------------------------------
    # Step 3
    # --------------------------------------------------------

    print(
        "\n[2] Starting docking..."
    )

    docking_results = perform_screening(
        ligands
    )

    # --------------------------------------------------------
    # Step 4
    # --------------------------------------------------------

    docking_results.to_csv(
        "data/results/docking_results.csv",
        index=False
    )

    print(
        "\nDocking results saved."
    )

    # --------------------------------------------------------
    # Step 5
    # --------------------------------------------------------

    matrix = build_binding_matrix(
        docking_results
    )

    matrix.to_csv(
        "data/results/"
        "conformation_binding_matrix.csv"
    )

    print(
        "\nConformation binding matrix:"
    )

    print(matrix.head())

    # --------------------------------------------------------
    # Step 6
    # --------------------------------------------------------

    analysis = (
        calculate_conformation_selectivity(
            matrix
        )
    )

    analysis.to_csv(
        "data/results/"
        "conformation_analysis.csv"
    )

    print(
        "\nAnalysis completed."
    )

    print(
        "\nTop compounds:"
    )

    print(
        analysis.sort_values(
            "best_score"
        ).head(20)
    )


if __name__ == "__main__":

    main()