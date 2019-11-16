require 'rest-client'
require 'net/http'
require 'bio'
require 'ruby-progressbar'
require './A3functions/functions'

# check number of input files and if name are duplicated
if ARGV.uniq.length != 1
  abort "\nError: wrong number of arguments.\n\nPlease use the command: ruby main.rb gene_codes_list.txt"
end

puts "\nWelcome to Assingment 3!"

puts "\nChecking files..."
agi_file = ARGV[0]
puts agi_file
fnameout = ARGV[1]
puts 'Succes!'

puts "\nFetching gene data from ENSEMBL..."
records = Functions.load_from_file(agi_file)

puts "\nWhich motif are you looking for?"
motif = $stdin.gets.chomp

puts "\nChecking genes that have the motif..."
genes_with_motif = Functions.genes_with_motif(records, motif)
puts "\nScanning exons..."
motif_coord_in_exon = Functions.scan_exons(genes_with_motif, motif)
puts "\nAdding features..."
Functions.add_motif_feature(genes_with_motif, motif_coord_in_exon, motif)

puts "\nChecking genes that did not have the motif..."
Functions.genes_without_motif(records, motif)
begin
  File.readlines("./output_files/genes_without_#{motif}.txt").each do |line|
    puts line
  end
rescue Errno::ENOENT
  puts 'Genes without motif were not found'
end

puts "\nCreating gff with gene coordinates..."
Functions.create_gff_genes(genes_with_motif, motif)

puts "\nCreating gff with genome coordinates..."
Functions.create_gff_genome(genes_with_motif, motif)