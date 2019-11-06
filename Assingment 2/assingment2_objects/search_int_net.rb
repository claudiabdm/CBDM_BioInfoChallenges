module Assingment2Objects
  ##
  # This class represents functional annotations.
  class AnnotationObject
    attr_accessor :database
    attr_accessor :entry_id
    attr_accessor :field
    attr_accessor :annotation

    ##
    # Creates a new annotation for a gene using the Togo REST API.
    def initialize(params = {})
      @database = params.fetch(:database, [])
      @entry_id = params.fetch(:entry_id, [])
      @field = params.fetch(:field, [])
      @annotation = load_annotation # automatically loads when new instances
    end

    ##
    # Fetch responses from the TOGO API using the values from the instance.
    def load_annotation
      begin
        address = "http://togows.org/entry/#{@database}/#{@entry_id}/#{@field}.json"
        response = RestClient::Request.execute(method: :get, url: address)
        @annotation = JSON.parse(response.body)[0] unless response.empty?
        return @annotation
      rescue RestClient::BadRequest => e
        abort e.response # prints the error from the Togo API
      end
      abort "\nError: Annotation not found" if @annotation.empty?
    end

    ##
    # Transform go annotation to hash[go_id] = go_term
    def go_to_hash(domain)
      @terms = {}
      unless annotation['GO'].nil?
        annotation['GO'].each do |term|
          next unless term[1].match(/#{domain}:/)

          key = term[0].split(':')[1]
          values = term[1]
          @terms[key] = values
        end
        annotation['GO'] = @terms
      end
    end
  end

  ##
  # Represent interactors that belong to gene list
  # and participates in PPI networks.
  class Interactor
    attr_accessor :agi_code
    attr_accessor :intact
    attr_accessor :kegg
    attr_accessor :go
    attr_accessor :list_interactors
    @@list_interactors = []

    ##
    # Creates a new interactor
    # instances are stored in a list automatically.
    def initialize(params = {})
      @agi_code = params.fetch(:agi_code, 'Not found')
      @intact = params.fetch(:intact, [])
      @kegg = params.fetch(:kegg, 'KEGG not found')
      @go = params.fetch(:go, 'GO not found')
      @@list_interactors << self
    end

    ##
    # Search the object that correponds to
    # that gene/interactor id (agi code)
    def self.get_interactor_info(id)
      id = String(id)
      interactor = @@list_interactors.find { |obj| obj.agi_code.casecmp(id).zero? }
      puts 'Error: gene id not found in database, please try again.' if interactor.nil?
      interactor
    end
  end

  ##
  # Represent interaction networks predicted from PPIs
  # in which genes from list are members of.
  class InteractionNetwork
    attr_accessor :network_num
    attr_accessor :predicted_net
    attr_accessor :members
    attr_accessor :go_common
    attr_accessor :kegg_all
    attr_accessor :list_possible_networks
    @@list_possible_networks = []

    ##
    # Creates new interaction newtwork that can be direct, 
    # gene from list interacts with gene from list and so on.
    # Or indirect, a middle interactor between genes from list
    def initialize(params = {})
      @network_num = params.fetch(:network_num, 0)
      @predicted_net = params.fetch(:predicted_net, 'Not found')
      @members = params.fetch(:members, 'Not found')
      @go_common = params.fetch(:go_common, 'Not found')
      @kegg_all = params.fetch(:kegg_all, 'Not found')
      @@list_possible_networks << self
    end

    ##
    # Search genes in IntAct from a gene codes list.
    # If response, genes are added to a hash that stores
    # the IntAct information (value) for the gene (key).
    def self.search_interactors(agi_file)
      @genes_intact = {} # stores genes that were found in intAct
      puts "\nSearching interactors from list in IntAct...\n"
      file_length = File.readlines(agi_file).length
      progressbar = ProgressBar.create(format: "%a %b\u{15E7}%i %p%% %t",
                                       progress_mark: ' ',
                                       remainder_mark: "\u{FF65}",
                                       starting_at: 0,
                                       total: file_length)
      # Searches every gene from the list in intAct
      File.open(agi_file).each do |line|
        address = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{line.strip}/?format=tab25"
        response = RestClient::Request.execute(method: :get, url: address)
        1.times { progressbar.increment }
        next if response.empty?

        # Store PPIs info for each gene from list
        @genes_intact[line.strip] = response.body.split("\n")
      end
      puts "\n#{@genes_intact.length} out of #{file_length} genes were found in IntAct."
      @genes_intact # return the list of genes found in IntAct
    end

    ##
    # Search PPIs from a hash whose keys are
    # the genes and values are the IntAct
    # information.
    # If filters are passed then PPIs are added
    # to another hash that contains all the genes
    # that genes from list interacts with.
    # Genes are filter by checking they belong to same species
    # and by the MIscore cutoff https://academic.oup.com/database/article/doi/10.1093/database/bau131/2433131
    def self.search_ppi(genes_intact, cutoff = 0.485, taxid = 'taxid:3702(arath)')
      puts "\nSearching protein-protein interactions (PPIs)..."
      @list_genes = genes_intact.keys
      @interacts_with = Hash.new { |h, k| h[k] = [] }
      count = 0
      # Filters
      species = /taxid:\d+\(\w+\)/
      agi_code = /(A[Tt][\d\w][Gg]\d\d\d\d\d)/
      miscore = /intact-miscore:(0\.\d+)/
      genes_intact.each do |_gene, intact|
        intact.each do |interaction|
          count += 1
          gene_a, gene_b = interaction.scan(agi_code)
          next if gene_a.nil? || gene_b.nil?

          species_a, species_b = interaction.scan(species)
          score = interaction.match(miscore).captures[0].to_f
          # PPIs added if species are the same and score is higher than cutoff
          next unless species_a == taxid && species_b == taxid && (score >= cutoff)

          # At least one gene is on the list
          if @list_genes.any? { |gene| gene.casecmp(gene_a[0]).zero? } || @list_genes.any? { |gene| gene.casecmp(gene_b[0]).zero? }
            @interacts_with[gene_a[0]] << gene_b[0]
            @interacts_with[gene_b[0]] << gene_a[0]
          end
        end
      end
      puts "Number of PPIs found: #{count}"
      puts "Number of genes that interacts with genes from list: #{@interacts_with.length} "
      @interacts_with
    end

    ##
    # Annotates with kegg and go functional annotation
    # the genes that belong to the list of interactors that
    # are part of some predicted network.
    def self.annot_interactors(final_interactors)
      puts "\nAnnotating interactors found, this might take a few minutes..."
      file_length = final_interactors.length
      progressbar = ProgressBar.create(format: "%a %b\u{15E7}%i %p%% %t",
                                       progress_mark: ' ',
                                       remainder_mark: "\u{FF65}",
                                       starting_at: 0,
                                       total: file_length)
      final_interactors.uniq.each do |gene|
        Interactor.new(agi_code: gene,
                       intact: @genes_intact[gene],
                       kegg: AnnotationObject.new(database: 'kegg-genes',
                                                  entry_id: "ath:#{gene}",
                                                  field: 'pathways').annotation,
                       go: AnnotationObject.new(database: 'ebi-uniprot',
                                                entry_id: gene.strip.to_s,
                                                field: 'dr').go_to_hash('P'))
        1.times { progressbar.increment }
      end
    end

    ##
    # Determines whether or not the interaction is direct or indirect
    # Direct (all gene are from gene list): A->B or A->B->C
    # Indirect (one gene is not from gene list): A -> B (not from list) -> C
    def self.network_type(interaction)
      if interaction.all? { |gene| @list_genes.any? { |gen| gen.casecmp(gene).zero? } }
        "#{interaction.join(' interacts with ')} (direct)"
      else
        "#{interaction.join(' interacts with ')} (indirect)"
      end
    end

    ##
    # Creates interaction networks where genes from the list interacts
    # with each other directly or through a middle interactor.
    # From a hash[gene_a] = gene_b determines interaction networks and
    # return a list of networks found where gene from list are members.
    def self.possible_networks(interactions)
      @interaction_networks = [] # store direct and indirect interactions
      interactions.keys.each do |gene_a|
        # the network starts with a gene that belongs to the list
        next if @list_genes.none? { |gene| gene.casecmp(gene_a).zero? }

        # check genes that interacts with first one
        # then the second one with a third one and stop
        interactions[gene_a].each do |gene_b|
          next if gene_a.casecmp(gene_b).zero? # genes must be different when direct connections

          interactions[gene_b].each do |gene_c|
            next if gene_c.casecmp(gene_a).zero? || gene_c.casecmp(gene_b).zero?

            # when third gene is not in the gene list, network is only composed by first and second (direct interaction)
            if @list_genes.none? { |gene| gene.casecmp(gene_c).zero? } && !gene_a.casecmp(gene_b).zero?
              arr = [gene_b, gene_a]
              next if @interaction_networks.any? { |network| network == arr } # try next if the reverse direction network already extist

              arr = [gene_a, gene_b]
              next if @interaction_networks.any? { |network| network == arr } # next if interaction already stored

              @interaction_networks << arr
              next
            end
            # networs three level deep
            arr = [gene_c, gene_b, gene_a]
            next if @interaction_networks.any? { |network| network == arr }

            arr = [gene_a, gene_b, gene_c]
            next if @interaction_networks.any? { |network| network == arr}

            @interaction_networks << arr
          end
        end
      end
      final_interactors = @interaction_networks.flatten.uniq # stores all genes used
      # Creates new interactor instances to annotate only the genes that at part of a network
      InteractionNetwork.annot_interactors(final_interactors)
      counter = 0
      @interaction_networks.each do |interaction|
        # Creates new object for each predicted network found and annotates the members
        InteractionNetwork.new(network_num: counter += 1,
                               predicted_net: InteractionNetwork.network_type(interaction),
                               members: interaction.map { |gene| Interactor.get_interactor_info(gene) },
                               go_common: Hash[*interaction.map { |gene| Interactor.get_interactor_info(gene) }.map(&:go).map(&:to_a).reduce(:&).flatten],
                               kegg_all: Hash[*interaction.map { |gene| Interactor.get_interactor_info(gene) }.map(&:kegg).map(&:to_a).flatten])
      end
      @@list_possible_networks # return a list with all InteractionNetwork objects created
    end
  end
end
