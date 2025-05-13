import tkinter as tk
from tkinter import PhotoImage, filedialog, messagebox

sequence = ""
query = ""
head_len = ""

def upload_file():
    global sequence
    file= filedialog.askopenfilename(
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa *.frn"), ("All files","*.*")]
    )
    if file:
        try:
            with open(file,"r") as file:
                global head_len
                lines = file.readlines()
                head = lines[0]
                head_len=len(head)
                print(f"{head}")
                sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
                dna_input.delete("1.0",tk.END)
                dna_input.insert("1.0",head+sequence)

        except Exception as e:
            messagebox.showerror("File Error", f"Error reading file: {e}")


def custom_fuzzy_search(seq, query_seq, match_val, mismatch_val, gap_penalty, top_hits):
    results = []
    seq = seq.upper()
    query_seq = query_seq.upper()

    for qlen in range(len(query_seq), 2, -1): #correction needed
        for qstart in range(0, len(query_seq) - qlen + 1):
            sub_query = query_seq[qstart:qstart + qlen]
            for i in range(len(seq) - qlen + 1):
                window = seq[i:i + qlen]
                score = 0
                for a, b in zip(sub_query, window):
                    if a == b:
                        score += match_val
                    else:
                        score += mismatch_val

                if score > 0:
                    results.append({
                        'qstart': qstart + 1,
                        'qend': qstart + qlen,
                        'query': sub_query,
                        'sstart': i + 1,
                        'send': i + qlen,
                        'subject': window,
                        'score': score
                    })
    sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)[:top_hits]
    return sorted_results

def search_func():
    global sequence
    query_seq = query_entry.get().strip().upper()
    head_sequence = dna_input.get("1.0", tk.END).strip().upper()
    sequence=head_sequence[head_len:]

    if not query_seq:
        messagebox.showerror("Warning", "Please enter a query pattern.")
        return
    if len(query_seq) > 30:
        messagebox.showerror("Warning", "Query pattern exceeds 30 characters.")
        return

    try:
        match = int(match_var.get())
        mismatch = int(mismatch_var.get())
        gap = int(gap_var.get())  # Placeholder, not used in current logic
        top_hits = int(report_var.get())
    except ValueError:
        messagebox.showerror("Error", "Invalid scoring parameters.")
        return

    results = custom_fuzzy_search(sequence, query_seq, match, mismatch, gap, top_hits)

    fuzz_output.delete("1.0", tk.END)
    if not results:
        fuzz_output.insert(tk.END, "No hits found.\n")
    else:
        fuzz_output.insert(tk.END, f"Fuzzy Search DNA results\nSearch results for {len(sequence)} residue sequence \"test sequence\" starting \"{sequence[:10].lower()}\"\nand {len(query_seq)} residue sequence \"query\" starting \"{query_seq[:10].lower()}\"\n\n")
        for res in results:
            fuzz_output.insert(tk.END, f">query from {res['qstart']} to {res['qend']}\n{res['query'].lower()}\n")
            fuzz_output.insert(tk.END, f">test sequence from {res['sstart']} to {res['send']}\n{res['subject'].lower()}\n")
            fuzz_output.insert(tk.END, f"Score: {res['score']}\n\n")

def reset_func():
    global sequence
    dna_input.delete("1.0",tk.END)
    query_entry.delete("0",tk.END)
    fuzz_output.delete("1.0",tk.END)
    match_var.set("1")
    mismatch_var.set("0")
    gap_var.set("-2")
    report_var.set("5")

app = tk.Tk()
app.title("BioFuzzAnalyzer")
app.resizable(False, False)
app.geometry("1000x600")
photo = PhotoImage(file="search (1).png")
app.iconphoto(False, photo)
app.configure(bg="white")

# --- Header ---
header = tk.Label(app,fg="white",bg="white")
header.pack(pady=(10,0))
image_label = tk.Label(header, image=photo,bg="white")
image_label.pack(side="left",padx=0, pady=0)
title_label = tk.Label(header, text="BioFuzzAnalyzer", font = ("Courier", 30, "bold"),bg="white",fg="#15BA77")
title_label.pack(side="right",padx=0, pady=0)

#Title
sub_title1 = tk.Label(app,font = ("Courier", 11), bg="white", text="BioFuzzAnalyzer enables efficient fuzzy DNA sequence searching,",fg="#37007F")
sub_title1.pack()
sub_title2 = tk.Label(app,font = ("Courier", 11), bg="white", text="allowing approximate matches to identify genetic variations and mutations accurately",fg="#37007F")
sub_title2.pack()

#Input
input_frame = tk.LabelFrame(text="Input Gene Sequence", bg="white",fg="#15BA77", font = ("Courier", 10, "bold"))
input_frame.pack(fill="x",padx=30,pady=10)
dna_input =tk.Text(input_frame,font=("courier",12),height=3)
dna_input.pack(fill="x",padx=5,pady=10)

#Option Frame
option_frame = tk.LabelFrame(text="Options", bg="white",fg="#15BA77", font = ("Courier", 10, "bold"))
option_frame.pack(fill="x",padx=30,pady=10)

#Option box 1
button1 = tk.Button(option_frame, text="Upload .FASTA file", bg="white", font = ("Courier", 10, "bold"), command=upload_file)
button1.pack(side="left", padx=10, pady= (5,10))
box1=tk.Label(option_frame, text="Query Pattern (Limit is 30 characters)",  font = ("Courier", 10, "bold"),bg="white",fg="#15BA77")
box1.pack(side="left", padx=10, pady=(5,10))
query_entry=tk.Entry(option_frame, width=45, font=("Courier", 11, "bold"),fg="#15BA77")
query_entry.pack(side="left", padx=10, pady=(5,10))

#Parameter frame
parameter_frame = tk.LabelFrame(text="Use the following parameters to specify how hits are scored, and to control which hits are returned.", bg="white",fg="#E55050", font = ("Courier", 10, "bold"))
parameter_frame.pack(fill="x",padx=30,pady=10)

#Parameter box 1
button2 = tk.Button(parameter_frame, text="Export", bg="white", font = ("Courier", 10, "bold"), width=6)
button2.pack(side="right", padx=(5,10), pady= (5,10))
button3 = tk.Button(parameter_frame, text="Reset", bg="white", font = ("Courier", 10, "bold"), width=6, command=reset_func)
button3.pack(side="right", padx=5, pady= (5,10))
button4 = tk.Button(parameter_frame, text="Search", bg="white", font = ("Courier", 10, "bold"), width=6, command=search_func)
button4.pack(side="right", padx=5, pady= (5,10))

#Dropdown variables and options
options = [str(i) for i in range (-5,5)]
match_var = tk.StringVar(value="1")
mismatch_var = tk.StringVar(value="0")
gap_var = tk.StringVar(value="-2")

report_options = [1]+[2]+[3]+[4]+[5]+[6]+[7]+[8]+[9]+[str(i) for i in range(10, 51, 10)]
report_var = tk.StringVar(value="5")

#Parameter box 2
box2=tk.Label(parameter_frame, text="Match Value",  font = ("Courier", 10, "bold"),bg="#15BA77",fg="white")
box2.pack(side="left", padx=(10,2), pady=(5,10))
match_entry=tk.OptionMenu(parameter_frame,match_var, *options)
match_entry.pack(side="left", padx=(4,5), pady=(5,10))

box3=tk.Label(parameter_frame, text="Mismatch Value",  font = ("Courier", 10, "bold"),bg="#15BA77",fg="white")
box3.pack(side="left", padx=2, pady=(5,10))
mismatch_entry=tk.OptionMenu(parameter_frame, mismatch_var, *options)
mismatch_entry.pack(side="left", padx=(4,5), pady=(5,10))

box4=tk.Label(parameter_frame, text="Gap Value",  font = ("Courier", 10, "bold"),bg="#15BA77",fg="white")
box4.pack(side="left", padx=2, pady=(5,10))
gap_entry=tk.OptionMenu(parameter_frame, gap_var, *options)
gap_entry.pack(side="left", padx=(4,5), pady=(5,10))

box5=tk.Label(parameter_frame, text="Report the top hits",  font = ("Courier", 10, "bold"),bg="#15BA77",fg="white")
box5.pack(side="left", padx=2, pady=(5,10))
report_entry=tk.OptionMenu(parameter_frame, report_var, *report_options)
report_entry.pack(side="left", padx=(4,5), pady=(5,10))

#Output Box
output_frame = tk.LabelFrame(text="Fuzz Analysis Report", bg="white",fg="#15BA77", font = ("Courier", 10, "bold"))
output_frame.pack(fill="x",padx=30,pady=10)
fuzz_output =tk.Text(output_frame,font=("courier",12),height=5)
fuzz_output.pack(fill="x",padx=5,pady=10)

#Copyright
copy_right = tk.Label(app,font = ("Courier", 8), bg="white", text="Developed by Jarif Ayman | © 2025 Bioinformatics Engineering, BAU",fg="#37007F")
copy_right.pack()


app.mainloop()