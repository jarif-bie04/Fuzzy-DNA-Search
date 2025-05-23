import tkinter as tk
from tkinter import PhotoImage, filedialog, messagebox
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from PIL import Image
from reportlab.lib.utils import ImageReader
from reportlab.lib import colors
from datetime import datetime
from os.path import basename

sequence = ""
query = ""
head_len = 0
fasta_filename = ""
fasta_identifier = ""

def upload_file():
    global sequence, fasta_filename, fasta_identifier
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
                fasta_filename = basename(file.name)

                if head.startswith(">"):
                    fasta_identifier = head[1:].split()[0]
                else:
                    fasta_identifier = "Unknown"

                sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
                dna_input.delete("1.0",tk.END)
                dna_input.insert("1.0",head+sequence)

        except Exception as e:
            messagebox.showerror("File Error", f"Error reading file: {e}")

def custom_fuzzy_search(seq, query_seq, match_val, mismatch_val, gap_penalty, top_hits):
    results = []
    seq = seq.upper()
    query_seq = query_seq.upper()

    for qlen in range(len(query_seq), 0, -1):
        for qstart in range(len(query_seq) - qlen + 1):
            sub_query = query_seq[qstart:qstart + qlen]

            for i in range(len(seq) - qlen + 1):
                subject = seq[i:i + qlen]

                m, n = len(sub_query), len(subject)
                dp = [[0] * (n + 1) for _ in range(m + 1)]

                for x in range(1, m + 1):
                    dp[x][0] = x * gap_penalty
                for y in range(1, n + 1):
                    dp[0][y] = y * gap_penalty

                for x in range(1, m + 1):
                    for y in range(1, n + 1):
                        match = dp[x - 1][y - 1] + (match_val if sub_query[x - 1] == subject[y - 1] else mismatch_val)
                        delete = dp[x - 1][y] + gap_penalty
                        insert = dp[x][y - 1] + gap_penalty
                        dp[x][y] = max(match, delete, insert)

                score = dp[m][n]
                results.append({
                    'qstart': qstart + 1,
                    'qend': qstart + qlen,
                    'query': sub_query,
                    'sstart': i + 1,
                    'send': i + qlen,
                    'subject': subject,
                    'score': score
                })

    sorted_results = sorted(results, key=lambda a: a['score'], reverse=True)[:top_hits]
    return sorted_results

def search_func():
    global sequence
    query_seq = query_entry.get().strip().upper()
    head_sequence = dna_input.get("1.0", tk.END).strip().upper()
    try:
        head_len_int = int(head_len)
    except ValueError:
        head_len_int = 0

    sequence = head_sequence[head_len_int:]

    if not query_seq:
        messagebox.showerror("Warning", "Please enter a query pattern.")
        return
    if len(query_seq) > 30:
        messagebox.showerror("Warning", "Query pattern exceeds 30 characters.")
        return

    try:
        match = int(match_var.get())
        mismatch = int(mismatch_var.get())
        gap = int(gap_var.get())
        top_hits = int(report_var.get())
    except ValueError:
        messagebox.showerror("Error", "Invalid scoring parameters.")
        return

    results = custom_fuzzy_search(sequence, query_seq, match, mismatch, gap, top_hits)

    fuzz_output.delete("1.0", tk.END)
    if not results:
        fuzz_output.insert(tk.END, "No hits found.\n")
    else:
        fuzz_output.insert(tk.END, f"Fuzzy Search DNA results\nSearch results for {len(sequence)} residue sequence \"{fasta_identifier}\" starting \"{sequence[:10].upper()}\"\nand {len(query_seq)} residue sequence \"query\" starting \"{query_seq[:10].upper()}\"\n\n")
        for res in results:
            fuzz_output.insert(tk.END, f">query from {res['qstart']} to {res['qend']}\n{res['query'].upper()}\n")
            fuzz_output.insert(tk.END, f">{fasta_identifier} from {res['sstart']} to {res['send']}\n{res['subject'].upper()}\n")
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

def export_to_pdf():
    try:
        file_path = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF files", "*.pdf")])
        if not file_path:
            return

        c = canvas.Canvas(file_path, pagesize=A4)
        width, height = A4
        y = height - 50

        export_date = datetime.now().strftime("Exported on: %Y-%m-%d %H:%M")
        c.setFont("Courier", 8)
        c.drawRightString(width - 40, height - 30, export_date)

        logo_path = "search (1).png"
        logo_width, logo_height = 80, 80
        try:
            with Image.open(logo_path) as img:
                if img.mode in ("RGBA", "LA"):
                    bg = Image.new("RGB", img.size, (255, 255, 255))
                    bg.paste(img, mask=img.split()[3])
                    logo_reader = ImageReader(bg)
                else:
                    logo_reader = ImageReader(img)

                logo_x = (width - logo_width) / 2
                c.drawImage(logo_reader, logo_x, y - logo_height, width=logo_width, height=logo_height)

        except Exception as e:
            print("Logo issue:", e)

        c.setFont("Courier", 16)
        c.setFillColor(colors.HexColor("#15BA77"))
        c.drawCentredString(width / 2, y - logo_height - 10, "BioFuzzAnalyzer")
        c.setFillColor(colors.black)

        y -= (logo_height + 40)

        c.setFont("Courier-Bold", 10)
        c.drawString(50, y, "Input Sequence:")
        y -= 15

        c.setFont("Courier", 10)
        input_seq = dna_input.get("1.0", tk.END).strip().splitlines()

        if input_seq:

            c.setFont("Courier", 10)
            header_line = input_seq[0]
            for i in range(0, len(header_line), 80):
                if y <= 40:
                    c.showPage()
                    y = height - 50
                    c.setFont("Courier", 10)
                c.drawString(50, y, header_line[i:i + 80])
                y -= 12

            y -=20

        if y <= 60:
            c.showPage()
            y = height - 50
        c.setFont("Courier-Bold", 10)
        c.drawString(50, y, "Fuzzy Search Output:")
        y -= 12

        c.setFont("Courier", 10)
        output_text = fuzz_output.get("1.0", tk.END).strip()
        for line in output_text.splitlines():
            if y <= 40:
                c.showPage()
                y = height - 50
                c.setFont("Courier", 10)
            c.drawString(50, y, line)
            y -= 12

        # Footer
        c.setFont("Courier", 8)
        c.setFillColor(colors.HexColor("#15BA77"))
        c.drawCentredString(width / 2, 30, "Developed by Jarif Ayman | © 2025 Bioinformatics Engineering, BAU")
        c.setFillColor(colors.black)

        c.save()
        messagebox.showinfo("Success", "PDF exported successfully")

    except Exception as e:
        messagebox.showerror("Export Error", f"Failed to export PDF: {e}")


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
button2 = tk.Button(parameter_frame, text="Export", bg="white", font = ("Courier", 10, "bold"), width=6, command=export_to_pdf)
button2.pack(side="right", padx=(5,10), pady= (5,10))
button3 = tk.Button(parameter_frame, text="Reset", bg="white", font = ("Courier", 10, "bold"), width=6, command=reset_func)
button3.pack(side="right", padx=5, pady= (5,10))
button4 = tk.Button(parameter_frame, text="Search", bg="white", font = ("Courier", 10, "bold"), width=6, command=search_func)
button4.pack(side="right", padx=5, pady= (5,10))

#Dropdown variables and options
options = [str(i) for i in range (-5,6)]
options_gap = [str(i) for i in range (-5,1)]
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
gap_entry=tk.OptionMenu(parameter_frame, gap_var, *options_gap)
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
