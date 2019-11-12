##
# Module for all function used in Assingment 3
module Functions
  module_function

  ##
  # Fetch information about gene in ensembl database.
  # Returns a hash with gene as key and info as value.
  def load_from_file(file)
    genes = File.readlines(file).map(&:strip) # removes new line \n
    records = {}
    genes.each do |gene|
      address = URI("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}")
      response = Net::HTTP.get_response(address)
      entry = Bio::EMBL.new(response.body)
      records[gene] = entry
    end
    puts "Genes loaded: #{records.length}"
    records
  end

  ##
  # Counts genes that have the repeat motif in both strands.
  # Input: hash[gene] = entry and repeat motif
  def scan_repeats(records, repeat_motif)
    search = Bio::Sequence::NA.new(repeat_motif.to_s)
    re = Regexp.new(search.to_re, Regexp::IGNORECASE)
    count = 0
    records.each do |_gene, entry|
      next if entry.to_biosequence.match(re).nil? && entry.to_biosequence.reverse_complement.match(re).nil?

      count += 1
    end
    puts "Genes were #{repeat_motif} sequence were found: #{count}"
  end

  ## 
  # Add a feature for the repeat motif in the exon for both strands.
  # Input: hash[gene] = entry and repeat motif.
  def add_features (records, repeat_motif)
    search = Bio::Sequence::NA.new(repeat_motif)
    re = Regexp.new(search)
    records.each do |_gene, entry|
      bioseq = entry.to_biosequence
      # Coordinates for the positive(seq_coord) and negative(comp_coord) strands
      # How to find the position: https://stackoverflow.com/a/51901945/10871520
      seq_coord = bioseq.gsub(re).map { [Regexp.last_match.begin(0), Regexp.last_match.end(0)-1]}
      comp_coord = bioseq.reverse_complement.gsub(re).map { [Regexp.last_match.begin(0), Regexp.last_match.end(0)-1]}
      # Creates and adds features to the ebsembl sequence object
      seq_coord.each do |coord|
        f1 = Bio::Feature.new('repeat_region', "#{coord[0]}..#{coord[1]}")
        f1.append(Bio::Feature::Qualifier.new('rpt_type', 'dispersed'))
        f1.append(Bio::Feature::Qualifier.new('function', 'target for insertions'))
        f1.append(Bio::Feature::Qualifier.new('rpt_unit_seq', repeat_motif.to_s))
        f1.append(Bio::Feature::Qualifier.new('strand', '+'))
        bioseq.features << f1
      end
      comp_coord.each do |coord|
        f2 = Bio::Feature.new('repeat_region', "complement(#{coord[0]}..#{coord[1]})")
        f2.append(Bio::Feature::Qualifier.new('rpt_type', 'dispersed'))
        f2.append(Bio::Feature::Qualifier.new('function', 'target for insertions'))
        f2.append(Bio::Feature::Qualifier.new('rpt_unit_seq', repeat_motif.to_s))
        f2.append(Bio::Feature::Qualifier.new('strand', '-'))
        bioseq.features << f2
      end
    end
    puts "\nFeatures added"
  end
  ##
  # Creates gff using gene coordinates.
  def create_gff_gene(records, repeat_motif, qualifier)
    File.open("genes_#{repeat_motif}.gff", 'w+') do |f|
      f.puts '##gff-version 3'
      records.each do |gene, entry|
        seqid = gene
        source = '.'
        score = '.'
        phase = '.'
        count = 0
        entry.features.each do |feature|
          # Filter features to be added
          next unless feature.feature == 'repeat_region' && feature.assoc['rpt_unit_seq'] == repeat_motif.to_s

          type = "#{feature.assoc[qualifier.to_s]}_repeat"
          start, ends = feature.position.match(/(\d+)\.\.(\d+)/).captures.map(&:to_i) # feature coordinates
          strand = feature.assoc['strand'].strip
          attributes = "ID=#{gene}_#{count+=1}_#{repeat_motif};start=#{start};end=#{ends}"
          # Add records to GFF file
          f.puts "#{seqid}\t#{source}\t#{type}\t#{start}\t#{ends}\t#{score}\t#{strand}\t#{phase}\t#{attributes}"
        end
      end
    end
    puts "genes_#{repeat_motif}.gff created and saved in #{Dir.pwd}\n"
  end

  ##
  # Creates a report with genes that do not have the repeat motif.
  def genes_without_repeat(records, repeat_motif)
    search = Bio::Sequence::NA.new(repeat_motif.to_s)
    re = Regexp.new(search.to_re, Regexp::IGNORECASE)
    File.open("genes_without_#{repeat_motif}.txt", 'w+') do |f|
      f.puts "Genes were #{repeat_motif} sequences were not found:"
      records.each do |gene, entry|
        # next if CTTCTT sequence is not found in any strand of the exon
        next unless entry.to_biosequence.match(re).nil? && entry.to_biosequence.reverse_complement.match(re).nil?

        f.puts gene.to_s
      end
    end
  end

  ##
  # Creates a gff file with genome coordinates
  def create_gff_genome (records, repeat_motif, qualifier)
    File.open("genome_#{repeat_motif}.gff", 'w+') do |f|
      f.puts '##gff-version 3'
      records.each do |gene, entry|
        seqid = entry.entry_id.strip
        coord = entry.ac[0].split(':')[3..4].map(&:to_i) # gene coordinates
        source = '.'
        score = '.'
        phase = '.'
        attributes = '.'
        count = 0
        entry.features.each do |feature|
          next unless feature.feature == 'repeat_region' && feature.assoc['rpt_unit_seq'] == repeat_motif.to_s

          type = "#{feature.assoc[qualifier.to_s]}_repeat"
          start, ends = feature.position.match(/(\d+)\.\.(\d+)/).captures.map(&:to_i) # feature coordinates
          # Scales gene coordinates into genome coordinates
          if feature.assoc['strand'] == '+'
            coord_start = coord[0] + start
            coord_end = coord[0] + ends
          else
            coord_end = coord[1] - start
            coord_start = coord[1] - ends
          end
          # Add records to GFF file
          attributes = "ID=#{gene}_#{count+=1}_#{repeat_motif};start=#{start};end=#{ends}" # ID + coords in gene
          strand = feature.assoc['strand'].strip
          f.puts "#{seqid}\t#{source}\t#{type}\t#{coord_start}\t#{coord_end}\t#{score}\t#{strand}\t#{phase}\t#{attributes}"
        end
      end
    end
    puts "genome_#{repeat_motif}.gff created and saved in #{Dir.pwd}\n"
  end
end
