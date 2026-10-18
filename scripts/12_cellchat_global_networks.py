from _bootstrap import finish_figure, pd, plt, setup, write_table


def main():
    _, repo, out = setup("Global CellChat network summaries", "12_cellchat_global")
    major = repo.table(26, "ST26a")
    fine = repo.table(26, "ST26b")
    write_table(major, out / "major_cell_network_summary.tsv")
    write_table(fine, out / "fine_state_outgoing_activity.tsv")
    for column in ["Outgoing strength", "Incoming strength", "Outgoing interaction count"]:
        major[column] = pd.to_numeric(major[column], errors="coerce")
    fig, ax = plt.subplots()
    points = ax.scatter(major["Outgoing strength"], major["Incoming strength"], s=20 + major["Outgoing interaction count"] / 12, c=range(len(major)), cmap="tab20", edgecolor="white")
    for _, row in major.iterrows():
        ax.text(row["Outgoing strength"], row["Incoming strength"], str(row["Cell class"]), fontsize=6)
    ax.set(xlabel="Outgoing interaction strength", ylabel="Incoming interaction strength")
    finish_figure(fig, out / "Fig6a_major_cell_cellchat_strength")
    fine["Outgoing interaction count"] = pd.to_numeric(fine["Outgoing interaction count"], errors="coerce")
    fine = fine.nlargest(35, "Outgoing interaction count").sort_values("Outgoing interaction count")
    fig, ax = plt.subplots()
    ax.barh(fine["Cell subtype"].astype(str), fine["Outgoing interaction count"], color="#1F6D8C")
    ax.set(xlabel="Outgoing interaction count", ylabel="Fine cell state")
    finish_figure(fig, out / "Fig6b_fine_state_cellchat_strength", 7, 8)


if __name__ == "__main__":
    main()
