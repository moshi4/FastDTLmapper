from pathlib import Path

import pytest
from fastdtlmapper.goea.goea import GOEA


def test_goea_default_param_plot(
    goea_target_gene_list_file: Path,
    goea_all_gene_list_file: Path,
    goea_association_file: Path,
    go_basic_obo_file: Path,
    tmp_path: Path,
):
    """test GOEA default parameter plot (Only check run successfully)"""
    goea = GOEA(
        target_list_file=goea_target_gene_list_file,
        all_list_file=goea_all_gene_list_file,
        association_file=goea_association_file,
        obo_file=go_basic_obo_file,
    )
    output_prefix = tmp_path / "goea"
    goea_result_file_list = goea.run(output_prefix)
    for goea_result_file in goea_result_file_list:
        goea.plot(goea_result_file, goea_result_file.with_suffix(""))
    # Check goea result file exists
    for goea_result_file in goea_result_file_list:
        assert goea_result_file.is_file()


def test_goea_user_specified_param_plot(
    goea_target_gene_list_file: Path,
    goea_all_gene_list_file: Path,
    goea_association_file: Path,
    go_basic_obo_file: Path,
    tmp_path: Path,
):
    """test GOEA default parameter plot (Only check run successfully)"""
    goea = GOEA(
        target_list_file=goea_target_gene_list_file,
        all_list_file=goea_all_gene_list_file,
        association_file=goea_association_file,
        obo_file=go_basic_obo_file,
        pvalue_thr=0.3,
        plot_max_num=20,
        plot_format="pdf",
        plot_color="#1affdb",
        use_adjusted_pvalue=True,
    )
    output_prefix = tmp_path / "goea"
    goea_result_file_list = goea.run(output_prefix)
    for goea_result_file in goea_result_file_list:
        goea.plot(goea_result_file, goea_result_file.with_suffix(""))
    # Check goea result file exists
    for goea_result_file in goea_result_file_list:
        assert goea_result_file.is_file()


@pytest.mark.skip(reason="This test depends on network status.")
def test_download_obo(tmp_path: Path) -> Path:
    """test go-basic.obo file downaload"""
    go_basic_obo_file = tmp_path / "go-basic.obo"
    GOEA.download_obo(go_basic_obo_file)
    assert go_basic_obo_file.is_file()


def test_extract_significant_goea_result_pvalue(goea_result_file: Path):
    """test extract_significant_goea_result (use pvalue)"""
    pvalue_thr = 0.01
    use_adjusted_pvalue = False
    min_depth = 2
    over_df, under_df = GOEA.extract_significant_goea_result(
        goea_result_file, pvalue_thr, use_adjusted_pvalue, min_depth
    )
    # Over representation dataframe check
    assert len(over_df[over_df["enrichment"] != "e"]) == 0
    assert len(over_df[over_df["p_uncorrected"] >= pvalue_thr]) == 0
    assert len(over_df[over_df["depth"] < min_depth]) == 0
    # Under representation dataframe check
    assert len(under_df[under_df["enrichment"] != "p"]) == 0
    assert len(under_df[under_df["p_uncorrected"] >= pvalue_thr]) == 0
    assert len(under_df[under_df["depth"] < min_depth]) == 0


def test_extract_significant_goea_result_qvalue(goea_result_file: Path):
    """test extract_significant_goea_result (use qvalue)"""
    pvalue_thr = 0.1
    use_adjusted_pvalue = True
    min_depth = 3
    over_df, under_df = GOEA.extract_significant_goea_result(
        goea_result_file, pvalue_thr, use_adjusted_pvalue, min_depth
    )
    # Over representation dataframe check
    assert len(over_df[over_df["enrichment"] != "e"]) == 0
    assert len(over_df[over_df["p_fdr_bh"] >= pvalue_thr]) == 0
    assert len(over_df[over_df["depth"] < min_depth]) == 0
    # Under representation dataframe check
    assert len(under_df[under_df["enrichment"] != "p"]) == 0
    assert len(under_df[under_df["p_fdr_bh"] >= pvalue_thr]) == 0
    assert len(under_df[under_df["depth"] < min_depth]) == 0


def test_format_significant_goea_dataframe(goea_result_file: Path):
    """test format_significant_goea_dataframe"""
    over_df, under_df = GOEA.extract_significant_goea_result(
        goea_result_file, pvalue_thr=0.05, use_adjusted_pvalue=False
    )
    node_id, gain_or_loss = "N00X", "gain"
    over_df = GOEA.format_significant_goea_dataframe(over_df, node_id, gain_or_loss)
    # Check user specified column(NODE_ID, GAIN/LOSS) add
    assert over_df["NODE_ID"].unique()[0] == node_id
    assert over_df["GAIN/LOSS"].unique()[0] == gain_or_loss
    # Check OVER/UNDER column value change
    assert over_df["OVER/UNDER"].unique()[0] == "over"
    # Check format columns contents and order
    assert list(over_df.columns) == [
        "NODE_ID",
        "GAIN/LOSS",
        "GO_CATEGORY",
        "OVER/UNDER",
        "GO",
        "GO_NAME",
        "RATIO_IN_STUDY",
        "RATIO_IN_POP",
        "PVALUE",
        "DEPTH",
        "STUDY_COUNT",
        "BH_ADJUSTED_PVALUE",
        "STUDY_ITEMS",
    ]
