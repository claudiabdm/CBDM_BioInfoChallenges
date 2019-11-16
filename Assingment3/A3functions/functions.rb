##
# Module for all function used in Assingment 3
module Functions
  module_function

  ##
  # Fetch information about gene in ensembl database.
  # Returns a hash with gene as key and info as value.
  def load_from_file(file)
    genes = File.readlines(file).map(&:strip) # list gene and removes \n
    records = {}
    progressbar = ProgressBar.create(format: "%a %b\u{15E7}%i %p%% %t",
                                     progress_mark: ' ',
                                     remainder_mark: "\u{FF65}",
                                     starting_at: 0,
                                     total: genes.length)
    genes.each do |gene|
      address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}")
      response = Net::HTTP.get_response(address)
      1.times { progressbar.increment }
      entry = Bio::EMBL.new(response.body)
      records[gene] = entry
    end
    puts "Genes loaded: #{records.length}"
    records
  end

  ##
  # Deletes genes where motif was not found
  # Returns a new records hash without those genes
  def genes_with_motif(records, motif)
    re = Regexp.new(Bio::Sequence::NA.new(motif.to_s).to_re)
    genes_with_motif = records.dup
    genes_with_motif.each do |gene, entry|
      # First, find out the gene strand and then, check that the motif is in the sequence.
      if entry.features.find { |feat| feat.feature == 'gene' }.position.match('complement').nil? # positive strand
        genes_with_motif.delete(gene) if entry.to_biosequence.match(re).nil?
      else # negative strand
        genes_with_motif.delete(gene) if entry.to_biosequence.reverse_complement.match(re).nil?
      end
    end
    puts "Number of genes where #{motif} sequences were found: #{genes_with_motif.length}"
    genes_with_motif
  end

  ##
  # Creates a report with genes that do not have the repeat motif.
  def genes_without_motif(records, motif)
    re = Regexp.new(Bio::Sequence::NA.new(motif.to_s).to_re)
    File.open("./output_files/genes_without_#{motif}.txt", 'w+') do |f|
      f.puts "Genes where #{motif} sequences were not found:"
      records.each do |gene, entry|
        # next if CTTCTT sequence is found in any strand of the exon
        next unless entry.to_biosequence.match(re).nil? && entry.to_biosequence.reverse_complement.match(re).nil?

        f.puts gene.to_s
      end
    end
    puts "genes_without_#{motif}.gff created and saved in #{Dir.pwd}/output_files\n"
  end

  ##
  # Get coordinates where motifs occur within the exons of genes
  # Returns a list with those coordinates
  def motif_coordinates(motif, gene_length, exon_position, strand, sequence)
    re = Regexp.new(Bio::Sequence::NA.new(motif.to_s).to_re)
    coordinates = []
    # Coordinates for the positive(seq_coord) and negative(comp_coord) strands
    # How to find the position: https://stackoverflow.com/a/51901945/10871520
    if strand.nil?
      motif_coord = sequence.gsub(re).map { Range.new(*[Regexp.last_match.begin(0), Regexp.last_match.end(0)-1])}
      [*motif_coord].each do |coord|
        coordinates << coord if (coord.first <= exon_position.last) && (exon_position.first <= coord.last)
      end
    else
      motif_coord = sequence.reverse_complement.gsub(re).map { Range.new(*[gene_length - Regexp.last_match.begin(0), gene_length - Regexp.last_match.end(0) +1 ])}
      [*motif_coord].each do |coord|
        coordinates << "complement(#{coord})" if (coord.first <= exon_position.last) && (exon_position.first <= coord.last)
      end
    end
    coordinates
  end

  ##
  # Scan motifs within exons of a gene
  # Return a hash where the key is the exon id and values are the motif coordinates
  def scan_exons(genes_with_motif, motif)
    re = Regexp.new(Bio::Sequence::NA.new(motif.to_s).to_re)
    motif_coord_in_exon = Hash.new { |h, k| h[k] = [] }
    genes_with_motif.each do |gene, entry|
      # Select the exon feature and check if motifs are in exons
      gene_position = entry.ac[0].split(':')[3..4].map(&:to_i) # gene coordinates
      gene_length = gene_position[1] - gene_position[0] # this is use to get the correct coordinates for complement sequence
      bioseq = entry.to_biosequence
      entry.features.each do |feature|
        next unless feature.feature == 'exon'

        exon_id = feature.assoc['note'].split('=')[1]
        next if exon_id.match(/#{gene}/i).nil?

        exon_position = Range.new(*feature.position.match(/(\d+)\.\.(\d+)/).captures.join('..').split('..').map(&:to_i)) # exon coordinates
        strand = feature.position.scan(/(\w+)\(\d+..\d+\)/).flatten[0] # forward (nil) or reverse strand (complement)
        coord = motif_coordinates(motif, gene_length, exon_position, strand, bioseq)
        motif_coord_in_exon[exon_id] << coord unless coord.empty? # if motif were not found, coordinates are nil
      end
    end
    total_exons = motif_coord_in_exon.keys.uniq.length
    puts "Number of exons where #{motif} sequences were found: #{total_exons}"
    motif_coord_in_exon
  end

  ##
  # Creates and adds repeat_region feature to the ensemble sequence object
  # The feature position is the motif coordinates
  def add_motif_feature(genes_with_motif, motif_coord_in_exon, motif)
    genes_with_motif.each do |gene, entry|
      bioseq = entry.to_biosequence
      motif_coord_in_exon.each do |exon, coords|
        next if exon.match(/#{gene}/i).nil? # to avoid exons be added in wrong gene

        coords.flatten.each do |coord|
          f1 = Bio::Feature.new('repeat_region', coord.to_s)
          f1.append(Bio::Feature::Qualifier.new('exon_id', exon))
          f1.append(Bio::Feature::Qualifier.new('rpt_type', 'dispersed'))
          f1.append(Bio::Feature::Qualifier.new('function', 'target for insertions'))
          f1.append(Bio::Feature::Qualifier.new('rpt_unit_seq', motif.to_s))
          if coord.to_s.match('complement')
            f1.append(Bio::Feature::Qualifier.new('strand', '-'))
          else
            f1.append(Bio::Feature::Qualifier.new('strand', '+'))
          end
          bioseq.features << f1
        end
      end
    end
    puts 'Features added'
  end

  ##
  # Creates a gff file with genes coordinates
  def create_gff_genes (genes_with_motif, motif)
    File.open("./output_files/genes_#{motif}.gff", 'w+') do |f|
      f.puts '##gff-version 3'
      exons_added = []
      genes_with_motif.each do |gene, entry|
        seqid = gene
        source = '.'
        score = '.'
        phase = '.'
        entry.features.each do |feature|
          next unless feature.assoc['rpt_unit_seq'] == motif.to_s

          type = "#{feature.assoc['rpt_type']}_repeat"
          strand = feature.assoc['strand'].strip
          exon_id = feature.assoc['exon_id']
          count = exons_added.count(exon_id) + 1 
          start, ends = feature.position.match(/(\d+)\.\.(\d+)/).captures.map(&:to_i) # feature coordinates
          attributes = "ID=#{motif}_#{exon_id}_#{count}; Parent=#{exon_id}; Note=#{feature.position}"
          # Add records to GFF file
          f.puts "#{seqid}\t#{source}\t#{type}\t#{start}\t#{ends}\t#{score}\t#{strand}\t#{phase}\t#{attributes}"
          # Register exon
          exons_added << exon_id
        end
      end
    end
    puts "genes_#{motif}.gff created and saved in #{Dir.pwd}/output_files\n"
  end

  ##
  # Creates a gff file with genome coordinates
  def create_gff_genome (genes_with_motif, motif)
    File.open("./output_files/genome_#{motif}.gff", 'w+') do |f|
      f.puts '##gff-version 3'
      exons_added = []
      genes_with_motif.each do |gene, entry|
        seqid = entry.entry_id.strip
        coord = entry.ac[0].split(':')[3..4].map(&:to_i) # gene coordinates
        source = '.'
        score = '.'
        phase = '.'
        entry.features.each do |feature|
          next unless feature.assoc['rpt_unit_seq'] == motif.to_s

          type = "#{feature.assoc['rpt_type']}_repeat"
          exon_id = feature.assoc['exon_id']
          count = exons_added.count(exon_id) + 1 
          strand = feature.assoc['strand'].strip
          start, ends = feature.position.match(/(\d+)\.\.(\d+)/).captures.map(&:to_i) # feature coordinates
          attributes = "ID=#{motif}_#{exon_id}_#{count}; Parent=#{exon_id}; Note=#{feature.position}" # ID + exon id + coords in gene
          # Scales gene coordinates into genome coordinates
          if strand == '+'
            coord_start = coord[0] + start
            coord_end = coord[0] + ends
          else
            coord_start = coord[0] + ends
            coord_end = coord[0] + start
          end
          # Add records to GFF file
          f.puts "#{seqid}\t#{source}\t#{type}\t#{coord_start}\t#{coord_end}\t#{score}\t#{strand}\t#{phase}\t#{attributes}"
          # Register exon
          exons_added << exon_id
        end
      end
    end
    puts "genome_#{motif}.gff created and saved in #{Dir.pwd}/output_files\n"
  end
end
